"""Property-based tests for the karyotype parser and validator.

Hypothesis generates structurally valid and deliberately broken karyotypes and
checks invariants that must hold for every input, complementing the
example-based tests in ``test_parser.py``.

Generated cases are *not* written to ``fixtures/validity.json``: that corpus is
shared verbatim with the TypeScript runner in the iscn-authenticator monorepo,
so any change to it must land in both repos together. The fixture-driven test
at the bottom only reads the corpus.

Hypothesis is a test-only dependency (``pip install -e .[test]``); the module is
skipped when it is not installed so ``python -m unittest discover tests`` keeps
working with the stdlib alone.
"""

import contextlib
import json
import pathlib
import unittest

try:
    from hypothesis import HealthCheck, assume, example, given, settings
    from hypothesis import strategies as st
except ImportError as exc:  # pragma: no cover - exercised only without the extra
    raise unittest.SkipTest("hypothesis not installed; pip install -e .[test]") from exc

from iscn_authenticator.main import validate_karyotype
from iscn_authenticator.parser import KaryotypeParser, ParseError

FIXTURES_PATH = pathlib.Path(__file__).resolve().parent.parent / "fixtures" / "validity.json"

settings.register_profile(
    "iscn",
    max_examples=200,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
settings.load_profile("iscn")

PARSER = KaryotypeParser()

# ---------------------------------------------------------------------------
# Strategies: building blocks
# ---------------------------------------------------------------------------

autosomes = st.integers(min_value=1, max_value=22).map(str)
chromosomes = st.one_of(autosomes, st.sampled_from(["X", "Y"]))
arms = st.sampled_from(["p", "q"])


@st.composite
def bands(draw, arm=None):
    """A breakpoint such as ``q13``, ``p11.2`` or ``q21.31``."""
    arm = arm or draw(arms)
    region_band = draw(st.integers(min_value=1, max_value=36))
    text = f"{arm}{region_band}"
    if draw(st.booleans()):
        text += "." + str(draw(st.integers(min_value=1, max_value=33)))
    return text


@st.composite
def same_arm_pair(draw):
    arm = draw(arms)
    return draw(bands(arm)) + draw(bands(arm))


@st.composite
def different_arm_pair(draw):
    return draw(bands("p")) + draw(bands("q"))


# ---------------------------------------------------------------------------
# Strategies: valid abnormalities, each paired with its expected (type, chromosome)
# ---------------------------------------------------------------------------


@st.composite
def numerical(draw):
    sign = draw(st.sampled_from(["+", "-"]))
    chrom = draw(chromosomes)
    return f"{sign}{chrom}", sign, chrom


@st.composite
def deletion_or_duplication(draw):
    kind = draw(st.sampled_from(["del", "dup"]))
    chrom = draw(chromosomes)
    bp = draw(st.one_of(bands(), same_arm_pair()))
    return f"{kind}({chrom})({bp})", kind, chrom


@st.composite
def inversion(draw):
    chrom = draw(chromosomes)
    bp = draw(st.one_of(same_arm_pair(), different_arm_pair()))
    return f"inv({chrom})({bp})", "inv", chrom


@st.composite
def translocation(draw):
    c1, c2 = draw(st.lists(autosomes, min_size=2, max_size=2, unique=True))
    b1, b2 = draw(bands()), draw(bands())
    return f"t({c1};{c2})({b1};{b2})", "t", f"{c1};{c2}"


@st.composite
def ring(draw):
    chrom = draw(chromosomes)
    return f"r({chrom})({draw(different_arm_pair())})", "r", chrom


@st.composite
def isochromosome(draw):
    chrom = draw(chromosomes)
    return f"i({chrom})({draw(bands())})", "i", chrom


valid_abnormalities = st.one_of(
    numerical(),
    deletion_or_duplication(),
    inversion(),
    translocation(),
    ring(),
    isochromosome(),
)

# A normal base with a coherent count/sex pair (the only combinations the
# SEX_CHR_COHERENCE rule constrains when no abnormalities are listed).
normal_bases = st.sampled_from(["46,XX", "46,XY", "45,X", "47,XXY", "47,XXX", "47,XYY", "46,U"])


@st.composite
def valid_karyotypes(draw):
    """Return (karyotype, [(type, chromosome), ...])."""
    base = draw(normal_bases)
    abns = draw(st.lists(valid_abnormalities, max_size=3))
    text = ",".join([base] + [a[0] for a in abns])
    return text, [(a[1], a[2]) for a in abns]


# ---------------------------------------------------------------------------
# Properties: valid inputs
# ---------------------------------------------------------------------------


class TestGeneratedValidKaryotypes(unittest.TestCase):
    @given(valid_karyotypes())
    def test_generated_karyotypes_validate(self, case):
        text, _ = case
        result = validate_karyotype(text)
        self.assertTrue(result.valid, f"{text!r}: {result.errors}")
        self.assertEqual(result.errors, [])
        self.assertIsNotNone(result.parsed)

    @given(valid_karyotypes())
    def test_ast_reflects_components(self, case):
        text, expected = case
        ast = PARSER.parse(text)
        count, sex = text.split(",")[:2]
        self.assertEqual(ast.chromosome_count, int(count))
        self.assertEqual(ast.sex_chromosomes, sex)
        self.assertEqual([(a.type, a.chromosome) for a in ast.abnormalities], expected)
        self.assertEqual(",".join([count, sex] + [a.raw for a in ast.abnormalities]), text)

    @given(valid_karyotypes(), st.text(" \t\n", max_size=3), st.text(" \t\n", max_size=3))
    def test_surrounding_whitespace_is_ignored(self, case, left, right):
        text, _ = case
        self.assertEqual(PARSER.parse(left + text + right), PARSER.parse(text))

    @given(
        st.lists(
            st.tuples(valid_karyotypes(), st.integers(min_value=1, max_value=99)),
            min_size=2,
            max_size=4,
        )
    )
    def test_mosaic_cell_lines(self, lines):
        text = "/".join(f"{k[0]}[{n}]" for k, n in lines)
        ast = PARSER.parse(text)
        self.assertEqual(len(ast.cell_lines), len(lines))
        self.assertEqual([cl.count for cl in ast.cell_lines], [n for _, n in lines])
        self.assertTrue(validate_karyotype(text).valid, text)

    @given(valid_abnormalities, st.sampled_from(["mat", "pat", "dn"]), st.booleans())
    def test_modifiers_are_stripped(self, abn, suffix, uncertain):
        token, abn_type, chrom = abn
        raw = ("?" if uncertain else "") + token + suffix
        [parsed] = PARSER.parse(f"46,XX,{raw}").abnormalities
        self.assertEqual((parsed.type, parsed.chromosome), (abn_type, chrom))
        self.assertEqual(parsed.inheritance, suffix)
        self.assertEqual(parsed.uncertain, uncertain)
        self.assertEqual(parsed.raw, raw)


# ---------------------------------------------------------------------------
# Properties: invalid inputs
# ---------------------------------------------------------------------------


class TestGeneratedInvalidKaryotypes(unittest.TestCase):
    @given(st.text())
    @example("\u00b2,XX")  # superscript two: str.isdigit() is True but int() raises
    def test_arbitrary_text_never_crashes(self, text):
        """The parser only ever raises ParseError; validate_karyotype never raises."""
        with contextlib.suppress(ParseError):
            PARSER.parse(text)
        result = validate_karyotype(text)
        self.assertIsInstance(result.valid, bool)
        if result.valid:
            self.assertEqual(result.errors, [])
            self.assertIsNotNone(result.parsed)
        else:
            self.assertTrue(result.errors)

    @given(st.text(min_size=1).filter(lambda s: s.strip() and "~" not in s and "/" not in s and "," not in s))
    def test_non_numeric_count_rejected(self, count):
        assume(not count.strip().isdigit())
        self.assertFalse(validate_karyotype(f"{count},XX").valid)

    @given(st.sampled_from(["\u00b2", "\u0664\u0666", "\uff14\uff16"]))
    def test_non_ascii_digit_count_rejected(self, count):
        """Only ASCII digits are accepted, as in the TypeScript port's /^\\d+$/."""
        self.assertFalse(validate_karyotype(f"{count},XX").valid)

    @given(st.one_of(st.integers(min_value=0, max_value=22), st.integers(min_value=93, max_value=10_000)))
    def test_count_out_of_range_rejected(self, count):
        self.assertFalse(validate_karyotype(f"{count},XX").valid)

    @given(st.text(alphabet="XYUABZxy0123456789 ", min_size=1, max_size=5))
    def test_bad_sex_chromosomes_rejected(self, sex):
        assume(set(sex.strip()) - set("XYU") or not sex.strip())
        self.assertFalse(validate_karyotype(f"46,{sex}").valid)

    @given(st.text(alphabet="Y", min_size=1, max_size=3))
    def test_sex_without_x_rejected(self, sex):
        self.assertFalse(validate_karyotype(f"46,{sex}").valid)

    @given(valid_karyotypes())
    def test_missing_comma_rejected(self, case):
        text, _ = case
        self.assertFalse(validate_karyotype(text.replace(",", "")).valid)

    @given(st.sampled_from(["del", "dup", "inv", "r", "trp", "qdp"]), chromosomes)
    def test_breakpoint_required(self, kind, chrom):
        self.assertFalse(validate_karyotype(f"46,XX,{kind}({chrom})").valid)

    @given(st.sampled_from(["del", "dup", "inv", "t"]), chromosomes, bands())
    def test_malformed_breakpoint_rejected(self, kind, chrom, band):
        token = f"t({chrom};1)({band}x;q1)" if kind == "t" else f"{kind}({chrom})({band}x)"
        self.assertFalse(validate_karyotype(f"46,XX,{token}").valid)

    @given(st.sampled_from(["del", "dup"]), chromosomes, bands("p"), bands("q"))
    def test_interstitial_across_arms_rejected(self, kind, chrom, p_band, q_band):
        self.assertFalse(validate_karyotype(f"46,XX,{kind}({chrom})({p_band}{q_band})").valid)

    @given(st.lists(autosomes, min_size=2, max_size=3, unique=True), bands())
    def test_translocation_breakpoint_count_must_match(self, chroms, band):
        self.assertFalse(validate_karyotype(f"46,XX,t({';'.join(chroms)})({band})").valid)


# ---------------------------------------------------------------------------
# Shared corpus: read-only, kept byte-identical with the TypeScript runner
# ---------------------------------------------------------------------------


class TestFixtureInvariants(unittest.TestCase):
    """Apply the whitespace invariant to every case in the shared corpus."""

    @classmethod
    def setUpClass(cls):
        with FIXTURES_PATH.open() as f:
            cls.fixtures = json.load(f)

    @given(st.data())
    def test_fixture_validity_survives_padding(self, data):
        cases = [(c["input"], True) for c in self.fixtures["valid"]]
        cases += [(c["input"], False) for c in self.fixtures["invalid"]]
        text, expected = data.draw(st.sampled_from(cases))
        pad = data.draw(st.text(" \t", max_size=3))
        self.assertIs(validate_karyotype(pad + text + pad).valid, expected)


if __name__ == "__main__":
    unittest.main()
