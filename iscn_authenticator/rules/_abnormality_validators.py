# iscn_authenticator/rules/_abnormality_validators.py
"""Validator implementations for abnormality rules.

This module is private. Tests and external callers should import
``Rule`` instances from ``iscn_authenticator.rules.abnormality`` instead.
"""

from collections.abc import Callable

from iscn_authenticator.models import Abnormality, KaryotypeAST

Validator = Callable[[KaryotypeAST, Abnormality], list[str]]

VALID_CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y"}

_COUNT_WORD = {0: "no", 1: "one", 2: "two", 3: "three"}


def _bp_phrase(counts: tuple[int, ...]) -> str:
    if counts == (0, 1):
        return "zero or one breakpoint"
    if counts == (1, 2):
        return "one or two breakpoints"
    if len(counts) == 1:
        n = counts[0]
        return f"{_COUNT_WORD[n]} breakpoint{'' if n == 1 else 's'}"
    return " or ".join(_COUNT_WORD[n] for n in counts) + " breakpoints"


def bp_count(
    abn_type: str,
    counts: int | tuple[int, ...],
    label: str,
    *,
    phrase: str | None = None,
) -> Validator:
    """Validator: type ``abn_type`` must have one of ``counts`` breakpoints."""
    allowed = (counts,) if isinstance(counts, int) else tuple(counts)
    phrase = phrase or _bp_phrase(allowed)

    def validate(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
        if abnormality.type != abn_type:
            return []
        bp = len(abnormality.breakpoints)
        if bp in allowed:
            return []
        return [f"{label} requires {phrase}, found {bp} in {abnormality.raw}"]

    return validate


def bp_matches_chromosomes(abn_type: str, label: str) -> Validator:
    """Validator: breakpoint count must equal semicolon-split chromosome count."""

    def validate(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
        if abnormality.type != abn_type:
            return []
        chr_count = len(abnormality.chromosome.split(";"))
        bp_count_ = len(abnormality.breakpoints)
        if chr_count == bp_count_:
            return []
        return [f"{label} has {chr_count} chromosomes but {bp_count_} breakpoints in {abnormality.raw}"]

    return validate


def numerical_chromosome(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Numerical (+/-) abnormalities must reference 1-22, X, or Y."""
    if abnormality.type not in ("+", "-"):
        return []
    if abnormality.chromosome in VALID_CHROMOSOMES:
        return []
    return [f"Invalid chromosome '{abnormality.chromosome}' in {abnormality.raw}. Must be 1-22, X, or Y"]


def breakpoint_arm(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Breakpoint arms must be ``p`` or ``q``."""
    if abnormality.type in ("+", "-", "unknown"):
        return []
    return [
        f"Invalid breakpoint arm '{bp.arm}' in {abnormality.raw}. Must be 'p' or 'q'"
        for bp in abnormality.breakpoints
        if bp.arm not in ("p", "q")
    ]


def _arms_match(abnormality: Abnormality, *, want_same: bool) -> tuple[str, str] | None:
    """Return (arm1, arm2) if the first two breakpoints violate the rule, else None."""
    arm1 = abnormality.breakpoints[0].arm
    arm2 = abnormality.breakpoints[1].arm
    same = arm1 == arm2
    if want_same and not same:
        return arm1, arm2
    if not want_same and same:
        return arm1, arm2
    return None


def _del_dup_validator(abn_type: str, label: str, arm_label: str) -> Validator:
    """Deletion/duplication: 1 or 2 breakpoints; if 2, both must be on the same arm."""
    count_check = bp_count(abn_type, (1, 2), label)

    def validate(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
        errors = count_check(ast, abnormality)
        if errors or abnormality.type != abn_type:
            return errors
        if len(abnormality.breakpoints) == 2:
            mismatch = _arms_match(abnormality, want_same=True)
            if mismatch:
                a1, a2 = mismatch
                return [f"{arm_label} breakpoints must be on same arm, found {a1} and {a2} in {abnormality.raw}"]
        return []

    return validate


def _two_bp_arm_validator(
    abn_type: str,
    label: str,
    *,
    want_same: bool,
) -> Validator:
    """Exactly two breakpoints on the same or different arms (depending on ``want_same``)."""
    count_check = bp_count(abn_type, 2, label)
    rule_phrase = "same arm" if want_same else "different arms"

    def validate(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
        errors = count_check(ast, abnormality)
        if errors or abnormality.type != abn_type:
            return errors
        mismatch = _arms_match(abnormality, want_same=want_same)
        if mismatch:
            a1, a2 = mismatch
            return [f"{label} breakpoints must be on {rule_phrase}, found {a1} and {a2} in {abnormality.raw}"]
        return []

    return validate


inversion = bp_count("inv", 2, "Inversion")
translocation = bp_matches_chromosomes("t", "Translocation")
deletion = _del_dup_validator("del", "Deletion", "Interstitial deletion")
duplication = _del_dup_validator("dup", "Duplication", "Duplication")
ring = _two_bp_arm_validator("r", "Ring chromosome", want_same=False)
isochromosome = bp_count("i", 1, "Isochromosome")
triplication = _two_bp_arm_validator("trp", "Triplication", want_same=True)
quadruplication = _two_bp_arm_validator("qdp", "Quadruplication", want_same=True)
dicentric = bp_matches_chromosomes("dic", "Dicentric")
isodicentric = bp_count("idic", 1, "Isodicentric")
robertsonian = bp_matches_chromosomes("rob", "Robertsonian translocation")
add_material = bp_count("add", 1, "Add")
fragile_site = bp_count("fra", 1, "Fragile site")
insertion = bp_count("ins", 3, "Insertion")
double_minutes = bp_count("dmin", 0, "Double minutes")
hsr = bp_count("hsr", (0, 1), "HSR")
marker = bp_count("mar", 0, "Marker chromosome")
pseudodicentric = bp_matches_chromosomes("psu dic", "Pseudodicentric")
acentric = bp_count("ace", (1, 2), "Acentric fragment", phrase="1-2 breakpoints")
telomeric_assoc = bp_matches_chromosomes("tas", "Telomeric association")
fission = bp_count("fis", 1, "Fission")
neocentromere = bp_count("neo", 1, "Neocentromere")
incomplete = bp_count("inc", 0, "Incomplete karyotype marker")
