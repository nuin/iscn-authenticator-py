# iscn_authenticator/_abnormality_parsers.py
"""Per-type abnormality parser methods, mixed into ``KaryotypeParser``."""

from iscn_authenticator._errors import ParseError
from iscn_authenticator._parser_patterns import (
    ACENTRIC_PATTERN,
    ADD_PATTERN,
    BREAKPOINT_PATTERN,
    DELETION_PATTERN,
    DERIVATIVE_PATTERN,
    DICENTRIC_PATTERN,
    DOUBLE_BREAKPOINT_PATTERN,
    DUPLICATION_PATTERN,
    FISSION_PATTERN,
    FRAGILE_SITE_PATTERN,
    HSR_LOCATION_PATTERN,
    INSERTION_PATTERN,
    INVERSION_PATTERN,
    ISOCHROMOSOME_LONG_PATTERN,
    ISOCHROMOSOME_SHORT_PATTERN,
    ISODICENTRIC_PATTERN,
    NEOCENTROMERE_PATTERN,
    PSEUDODICENTRIC_PATTERN,
    QUADRUPLICATION_PATTERN,
    RING_BREAKPOINT_PATTERN,
    RING_SIMPLE_PATTERN,
    ROBERTSONIAN_PATTERN,
    TELOMERIC_ASSOC_PATTERN,
    TRANSLOCATION_PATTERN,
    TRIPLE_BREAKPOINT_PATTERN,
    TRIPLICATION_PATTERN,
)
from iscn_authenticator.models import Abnormality, Breakpoint


def _abn(type_: str, chromosome: str, breakpoints: list[Breakpoint], raw: str) -> Abnormality:
    """Construct an Abnormality with default optional fields."""
    return Abnormality(
        type=type_,
        chromosome=chromosome,
        breakpoints=breakpoints,
        inheritance=None,
        uncertain=False,
        copy_count=None,
        raw=raw,
    )


class _AbnormalityParsers:
    """Mixin: per-type parsers used by ``KaryotypeParser``."""

    def _parse_breakpoint(self, bp_str: str) -> Breakpoint:
        """Parse a single breakpoint like 'q13' or 'p11.2'."""
        match = BREAKPOINT_PATTERN.match(bp_str)
        if not match:
            raise ParseError(f"Invalid breakpoint format: '{bp_str}'")

        arm, region_band, subband = match.group(1), match.group(2), match.group(3)
        if len(region_band) >= 2:
            region = int(region_band[0])
            band = int(region_band[1:])
        else:
            region = int(region_band)
            band = 0
        return Breakpoint(arm=arm, region=region, band=band, subband=subband, uncertain=False)

    def _parse_breakpoints(self, bp_str: str) -> list[Breakpoint]:
        """Parse one or two concatenated breakpoints like 'q21' or 'q21q31'."""
        double_bp = DOUBLE_BREAKPOINT_PATTERN.match(bp_str)
        if double_bp:
            return [
                self._parse_breakpoint(double_bp.group(1)),
                self._parse_breakpoint(double_bp.group(2)),
            ]
        return [self._parse_breakpoint(bp_str)]

    def _parse_multiple_breakpoints(self, bp_str: str) -> list[Breakpoint]:
        """Parse semicolon-separated breakpoints like 'p11;p11'."""
        return [self._parse_breakpoint(part.strip()) for part in bp_str.split(";")]

    def _parse_two_breakpoints(self, bp_str: str, kind: str) -> list[Breakpoint]:
        """Parse exactly two concatenated breakpoints; raise if not present."""
        double_bp = DOUBLE_BREAKPOINT_PATTERN.match(bp_str)
        if not double_bp:
            raise ParseError(f"{kind} requires two breakpoints: '{bp_str}'")
        return [
            self._parse_breakpoint(double_bp.group(1)),
            self._parse_breakpoint(double_bp.group(2)),
        ]

    @staticmethod
    def _match_or_raise(pattern, part: str, kind: str):
        match = pattern.match(part)
        if not match:
            raise ParseError(f"Invalid {kind} format: '{part}'")
        return match

    def _parse_deletion(self, part: str) -> Abnormality:
        m = self._match_or_raise(DELETION_PATTERN, part, "deletion")
        return _abn("del", m.group(1), self._parse_breakpoints(m.group(2)), part)

    def _parse_duplication(self, part: str) -> Abnormality:
        m = self._match_or_raise(DUPLICATION_PATTERN, part, "duplication")
        return _abn("dup", m.group(1), self._parse_breakpoints(m.group(2)), part)

    def _parse_inversion(self, part: str) -> Abnormality:
        m = self._match_or_raise(INVERSION_PATTERN, part, "inversion")
        return _abn("inv", m.group(1), self._parse_two_breakpoints(m.group(2), "Inversion"), part)

    def _parse_translocation(self, part: str) -> Abnormality:
        m = self._match_or_raise(TRANSLOCATION_PATTERN, part, "translocation")
        return _abn("t", m.group(1), self._parse_multiple_breakpoints(m.group(2)), part)

    def _parse_isochromosome(self, part: str) -> Abnormality:
        short = ISOCHROMOSOME_SHORT_PATTERN.match(part)
        if short:
            bp = Breakpoint(arm=short.group(2), region=1, band=0, subband=None, uncertain=False)
            return _abn("i", short.group(1), [bp], part)
        long_ = ISOCHROMOSOME_LONG_PATTERN.match(part)
        if long_:
            return _abn("i", long_.group(1), [self._parse_breakpoint(long_.group(2))], part)
        raise ParseError(f"Invalid isochromosome format: '{part}'")

    def _parse_ring(self, part: str) -> Abnormality:
        simple = RING_SIMPLE_PATTERN.match(part)
        if simple:
            return _abn("r", simple.group(1), [], part)
        bp_match = RING_BREAKPOINT_PATTERN.match(part)
        if bp_match:
            breakpoints: list[Breakpoint] = []
            double_bp = DOUBLE_BREAKPOINT_PATTERN.match(bp_match.group(2))
            if double_bp:
                breakpoints = [
                    self._parse_breakpoint(double_bp.group(1)),
                    self._parse_breakpoint(double_bp.group(2)),
                ]
            return _abn("r", bp_match.group(1), breakpoints, part)
        raise ParseError(f"Invalid ring chromosome format: '{part}'")

    def _parse_insertion(self, part: str) -> Abnormality:
        """Parse ins(5;2)(p14;q21q31) or ins(2)(p13q21q31)."""
        m = self._match_or_raise(INSERTION_PATTERN, part, "insertion")
        chromosomes_str, breakpoints_str = m.group(1), m.group(2)

        if ";" in breakpoints_str:
            site_str, segment_str = (s.strip() for s in breakpoints_str.split(";", 1))
            breakpoints = [self._parse_breakpoint(site_str)]
            double_bp = DOUBLE_BREAKPOINT_PATTERN.match(segment_str)
            if double_bp:
                breakpoints.append(self._parse_breakpoint(double_bp.group(1)))
                breakpoints.append(self._parse_breakpoint(double_bp.group(2)))
            else:
                breakpoints.append(self._parse_breakpoint(segment_str))
            return _abn("ins", chromosomes_str, breakpoints, part)

        triple = TRIPLE_BREAKPOINT_PATTERN.match(breakpoints_str)
        if not triple:
            raise ParseError(f"Invalid insertion breakpoints: '{breakpoints_str}'")
        return _abn(
            "ins",
            chromosomes_str,
            [self._parse_breakpoint(triple.group(i)) for i in (1, 2, 3)],
            part,
        )

    def _parse_add(self, part: str) -> Abnormality:
        m = self._match_or_raise(ADD_PATTERN, part, "additional material")
        return _abn("add", m.group(1), [self._parse_breakpoint(m.group(2))], part)

    def _parse_triplication(self, part: str) -> Abnormality:
        m = self._match_or_raise(TRIPLICATION_PATTERN, part, "triplication")
        bps = self._parse_two_breakpoints(m.group(2), "Triplication")
        return _abn("trp", m.group(1), bps, part)

    def _parse_dicentric(self, part: str) -> Abnormality:
        m = self._match_or_raise(DICENTRIC_PATTERN, part, "dicentric")
        return _abn("dic", m.group(1), self._parse_multiple_breakpoints(m.group(2)), part)

    def _parse_isodicentric(self, part: str) -> Abnormality:
        m = self._match_or_raise(ISODICENTRIC_PATTERN, part, "isodicentric")
        return _abn("idic", m.group(1), [self._parse_breakpoint(m.group(2))], part)

    def _parse_fragile_site(self, part: str) -> Abnormality:
        m = self._match_or_raise(FRAGILE_SITE_PATTERN, part, "fragile site")
        return _abn("fra", m.group(1), [self._parse_breakpoint(m.group(2))], part)

    def _parse_robertsonian(self, part: str) -> Abnormality:
        m = self._match_or_raise(ROBERTSONIAN_PATTERN, part, "Robertsonian translocation")
        return _abn("rob", m.group(1), self._parse_multiple_breakpoints(m.group(2)), part)

    def _parse_quadruplication(self, part: str) -> Abnormality:
        m = self._match_or_raise(QUADRUPLICATION_PATTERN, part, "quadruplication")
        bps = self._parse_two_breakpoints(m.group(2), "Quadruplication")
        return _abn("qdp", m.group(1), bps, part)

    def _parse_hsr_location(self, part: str) -> Abnormality:
        m = self._match_or_raise(HSR_LOCATION_PATTERN, part, "HSR")
        return _abn("hsr", m.group(1), [self._parse_breakpoint(m.group(2))], part)

    def _parse_pseudodicentric(self, part: str) -> Abnormality:
        m = self._match_or_raise(PSEUDODICENTRIC_PATTERN, part, "pseudodicentric")
        return _abn("psu dic", m.group(1), self._parse_multiple_breakpoints(m.group(2)), part)

    def _parse_derivative(self, part: str) -> Abnormality:
        m = self._match_or_raise(DERIVATIVE_PATTERN, part, "derivative")
        return _abn("der", m.group(1), [], part)

    def _parse_acentric(self, part: str) -> Abnormality:
        m = self._match_or_raise(ACENTRIC_PATTERN, part, "acentric")
        return _abn("ace", m.group(1), self._parse_breakpoints(m.group(2)), part)

    def _parse_telomeric_assoc(self, part: str) -> Abnormality:
        m = self._match_or_raise(TELOMERIC_ASSOC_PATTERN, part, "telomeric association")
        return _abn("tas", m.group(1), self._parse_multiple_breakpoints(m.group(2)), part)

    def _parse_fission(self, part: str) -> Abnormality:
        m = self._match_or_raise(FISSION_PATTERN, part, "fission")
        return _abn("fis", m.group(1), [self._parse_breakpoint(m.group(2))], part)

    def _parse_neocentromere(self, part: str) -> Abnormality:
        m = self._match_or_raise(NEOCENTROMERE_PATTERN, part, "neocentromere")
        return _abn("neo", m.group(1), [self._parse_breakpoint(m.group(2))], part)
