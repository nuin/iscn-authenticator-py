# iscn_authenticator/rules/abnormality.py
"""Validation rules for karyotype abnormalities."""

from iscn_authenticator.models import Abnormality, KaryotypeAST
from iscn_authenticator.rules.base import Rule

# Valid autosome numbers (1-22) and sex chromosomes (X, Y)
VALID_CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y"}


def _validate_numerical_chromosome(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate that numerical abnormalities reference valid chromosomes."""
    # Only apply to numerical abnormalities (+, -)
    if abnormality.type not in ("+", "-"):
        return []

    chromosome = abnormality.chromosome
    if chromosome not in VALID_CHROMOSOMES:
        return [f"Invalid chromosome '{chromosome}' in {abnormality.raw}. Must be 1-22, X, or Y"]
    return []


def _validate_breakpoint_arm(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate that breakpoint arms are p or q."""
    # Only apply to structural abnormalities with breakpoints
    if abnormality.type in ("+", "-", "unknown"):
        return []

    errors = []
    for bp in abnormality.breakpoints:
        if bp.arm not in ("p", "q"):
            errors.append(f"Invalid breakpoint arm '{bp.arm}' in {abnormality.raw}. Must be 'p' or 'q'")
    return errors


def _validate_inversion_two_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate that inversions have exactly two breakpoints."""
    if abnormality.type != "inv":
        return []

    bp_count = len(abnormality.breakpoints)
    if bp_count != 2:
        return [f"Inversion requires two breakpoints, found {bp_count} in {abnormality.raw}"]
    return []


def _validate_translocation_breakpoint_count(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate that translocation breakpoint count matches chromosome count."""
    if abnormality.type != "t":
        return []

    # Count chromosomes from the chromosome field (semicolon-separated)
    chromosomes = abnormality.chromosome.split(";")
    chr_count = len(chromosomes)
    bp_count = len(abnormality.breakpoints)

    if chr_count != bp_count:
        return [f"Translocation has {chr_count} chromosomes but {bp_count} breakpoints in {abnormality.raw}"]
    return []


def _validate_deletion_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate deletion breakpoints.

    - Terminal deletion: 1 breakpoint
    - Interstitial deletion: 2 breakpoints on same arm
    """
    if abnormality.type != "del":
        return []

    bp_count = len(abnormality.breakpoints)

    # Must have 1 or 2 breakpoints
    if bp_count not in (1, 2):
        return [f"Deletion requires one or two breakpoints, found {bp_count} in {abnormality.raw}"]

    # If 2 breakpoints, they must be on the same arm
    if bp_count == 2:
        arm1 = abnormality.breakpoints[0].arm
        arm2 = abnormality.breakpoints[1].arm
        if arm1 != arm2:
            return [
                f"Interstitial deletion breakpoints must be on same arm, found {arm1} and {arm2} in {abnormality.raw}"
            ]

    return []


def _validate_duplication_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate duplication breakpoints.

    - Tandem duplication: 1 breakpoint
    - Interstitial duplication: 2 breakpoints on same arm
    """
    if abnormality.type != "dup":
        return []

    bp_count = len(abnormality.breakpoints)

    # Must have 1 or 2 breakpoints
    if bp_count not in (1, 2):
        return [f"Duplication requires one or two breakpoints, found {bp_count} in {abnormality.raw}"]

    # If 2 breakpoints, they must be on the same arm
    if bp_count == 2:
        arm1 = abnormality.breakpoints[0].arm
        arm2 = abnormality.breakpoints[1].arm
        if arm1 != arm2:
            return [f"Duplication breakpoints must be on same arm, found {arm1} and {arm2} in {abnormality.raw}"]

    return []


def _validate_ring_chromosome_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate ring chromosome breakpoints.

    Ring chromosomes require exactly 2 breakpoints on different arms (p and q).
    """
    if abnormality.type != "r":
        return []

    bp_count = len(abnormality.breakpoints)

    # Must have exactly 2 breakpoints
    if bp_count != 2:
        return [f"Ring chromosome requires two breakpoints, found {bp_count} in {abnormality.raw}"]

    # Breakpoints must be on different arms
    arm1 = abnormality.breakpoints[0].arm
    arm2 = abnormality.breakpoints[1].arm
    if arm1 == arm2:
        return [f"Ring chromosome breakpoints must be on different arms, found {arm1} and {arm2} in {abnormality.raw}"]

    return []


def _validate_isochromosome_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate isochromosome breakpoints.

    Isochromosomes have exactly 1 breakpoint at the centromere region.
    """
    if abnormality.type != "i":
        return []

    bp_count = len(abnormality.breakpoints)

    # Must have exactly 1 breakpoint
    if bp_count != 1:
        return [f"Isochromosome requires one breakpoint, found {bp_count} in {abnormality.raw}"]

    return []


def _validate_triplication_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate triplication breakpoints.

    Triplications have exactly 2 breakpoints on the same arm.
    """
    if abnormality.type != "trp":
        return []

    bp_count = len(abnormality.breakpoints)

    # Must have exactly 2 breakpoints
    if bp_count != 2:
        return [f"Triplication requires two breakpoints, found {bp_count} in {abnormality.raw}"]

    # Breakpoints must be on the same arm
    arm1 = abnormality.breakpoints[0].arm
    arm2 = abnormality.breakpoints[1].arm
    if arm1 != arm2:
        return [f"Triplication breakpoints must be on same arm, found {arm1} and {arm2} in {abnormality.raw}"]

    return []


def _validate_quadruplication_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate quadruplication breakpoints.

    Quadruplications have exactly 2 breakpoints on the same arm.
    """
    if abnormality.type != "qdp":
        return []

    bp_count = len(abnormality.breakpoints)

    # Must have exactly 2 breakpoints
    if bp_count != 2:
        return [f"Quadruplication requires two breakpoints, found {bp_count} in {abnormality.raw}"]

    # Breakpoints must be on the same arm
    arm1 = abnormality.breakpoints[0].arm
    arm2 = abnormality.breakpoints[1].arm
    if arm1 != arm2:
        return [f"Quadruplication breakpoints must be on same arm, found {arm1} and {arm2} in {abnormality.raw}"]

    return []


def _validate_dicentric_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate dicentric chromosome breakpoints.

    Dicentric chromosomes have breakpoint count matching chromosome count.
    """
    if abnormality.type != "dic":
        return []

    # Count chromosomes from the chromosome field (semicolon-separated)
    chromosomes = abnormality.chromosome.split(";")
    chr_count = len(chromosomes)
    bp_count = len(abnormality.breakpoints)

    if chr_count != bp_count:
        return [f"Dicentric has {chr_count} chromosomes but {bp_count} breakpoints in {abnormality.raw}"]
    return []


def _validate_isodicentric_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate isodicentric chromosome breakpoints.

    Isodicentric chromosomes have exactly 1 breakpoint.
    """
    if abnormality.type != "idic":
        return []

    bp_count = len(abnormality.breakpoints)

    # Must have exactly 1 breakpoint
    if bp_count != 1:
        return [f"Isodicentric requires one breakpoint, found {bp_count} in {abnormality.raw}"]

    return []


def _validate_robertsonian_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate Robertsonian translocation breakpoints.

    Robertsonian translocations have breakpoint count matching chromosome count.
    """
    if abnormality.type != "rob":
        return []

    # Count chromosomes from the chromosome field (semicolon-separated)
    chromosomes = abnormality.chromosome.split(";")
    chr_count = len(chromosomes)
    bp_count = len(abnormality.breakpoints)

    if chr_count != bp_count:
        return [
            f"Robertsonian translocation has {chr_count} chromosomes but {bp_count} breakpoints in {abnormality.raw}"
        ]
    return []


def _validate_add_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate additional material (add) breakpoints.

    Add has exactly 1 breakpoint indicating where the unknown material is attached.
    """
    if abnormality.type != "add":
        return []

    bp_count = len(abnormality.breakpoints)

    if bp_count != 1:
        return [f"Add requires one breakpoint, found {bp_count} in {abnormality.raw}"]

    return []


def _validate_fra_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate fragile site (fra) breakpoints.

    Fragile site has exactly 1 breakpoint indicating the fragile location.
    """
    if abnormality.type != "fra":
        return []

    bp_count = len(abnormality.breakpoints)

    if bp_count != 1:
        return [f"Fragile site requires one breakpoint, found {bp_count} in {abnormality.raw}"]

    return []


def _validate_ins_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate insertion (ins) breakpoints.

    Insertions have 3 breakpoints: 1 for insertion site + 2 for the inserted segment.
    For insertions within the same chromosome, all 3 breakpoints are listed.
    For insertions between chromosomes, format is ins(A;B)(pX;pYpZ).
    """
    if abnormality.type != "ins":
        return []

    bp_count = len(abnormality.breakpoints)

    if bp_count != 3:
        return [f"Insertion requires three breakpoints, found {bp_count} in {abnormality.raw}"]

    return []


def _validate_dmin_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate double minutes (dmin) breakpoints.

    Double minutes have no breakpoints - they are extrachromosomal elements.
    """
    if abnormality.type != "dmin":
        return []

    bp_count = len(abnormality.breakpoints)

    if bp_count != 0:
        return [f"Double minutes should have no breakpoints, found {bp_count} in {abnormality.raw}"]

    return []


def _validate_hsr_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate homogeneously staining region (hsr) breakpoints.

    HSR can have 0 breakpoints (simple notation) or 1 breakpoint (with location).
    """
    if abnormality.type != "hsr":
        return []

    bp_count = len(abnormality.breakpoints)

    if bp_count > 1:
        return [f"HSR should have zero or one breakpoint, found {bp_count} in {abnormality.raw}"]

    return []


def _validate_mar_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate marker chromosome (mar) breakpoints.

    Marker chromosomes have no breakpoints - origin is unknown.
    """
    if abnormality.type != "mar":
        return []

    bp_count = len(abnormality.breakpoints)

    if bp_count != 0:
        return [f"Marker chromosome should have no breakpoints, found {bp_count} in {abnormality.raw}"]

    return []


def _validate_pseudodicentric_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate pseudodicentric (psu dic) breakpoints.

    Pseudodicentric chromosomes have breakpoint count matching chromosome count.
    """
    if abnormality.type != "psu dic":
        return []

    # Count chromosomes from the chromosome field (semicolon-separated)
    chromosomes = abnormality.chromosome.split(";")
    chr_count = len(chromosomes)
    bp_count = len(abnormality.breakpoints)

    if chr_count != bp_count:
        return [f"Pseudodicentric has {chr_count} chromosomes but {bp_count} breakpoints in {abnormality.raw}"]
    return []


def _validate_acentric_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate acentric fragment (ace) breakpoints.

    Acentric fragments have 1-2 breakpoints defining the fragment.
    """
    if abnormality.type != "ace":
        return []

    bp_count = len(abnormality.breakpoints)

    if bp_count not in (1, 2):
        return [f"Acentric fragment requires 1-2 breakpoints, found {bp_count} in {abnormality.raw}"]

    return []


def _validate_telomeric_association_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate telomeric association (tas) breakpoints.

    Telomeric associations have breakpoint count matching chromosome count.
    """
    if abnormality.type != "tas":
        return []

    # Count chromosomes from the chromosome field (semicolon-separated)
    chromosomes = abnormality.chromosome.split(";")
    chr_count = len(chromosomes)
    bp_count = len(abnormality.breakpoints)

    if chr_count != bp_count:
        return [f"Telomeric association has {chr_count} chromosomes but {bp_count} breakpoints in {abnormality.raw}"]
    return []


def _validate_fission_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate fission (fis) breakpoints.

    Fission has exactly 1 breakpoint at the centromere.
    """
    if abnormality.type != "fis":
        return []

    bp_count = len(abnormality.breakpoints)

    if bp_count != 1:
        return [f"Fission requires one breakpoint, found {bp_count} in {abnormality.raw}"]

    return []


def _validate_neocentromere_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate neocentromere (neo) breakpoints.

    Neocentromere has exactly 1 breakpoint indicating the new centromere location.
    """
    if abnormality.type != "neo":
        return []

    bp_count = len(abnormality.breakpoints)

    if bp_count != 1:
        return [f"Neocentromere requires one breakpoint, found {bp_count} in {abnormality.raw}"]

    return []


def _validate_incomplete_breakpoints(ast: KaryotypeAST, abnormality: Abnormality) -> list[str]:
    """Validate incomplete karyotype (inc) breakpoints.

    Incomplete marker has no breakpoints - it's just a notation marker.
    """
    if abnormality.type != "inc":
        return []

    bp_count = len(abnormality.breakpoints)

    if bp_count != 0:
        return [f"Incomplete karyotype marker should have no breakpoints, found {bp_count} in {abnormality.raw}"]

    return []


# Rule instances
numerical_chromosome_valid_rule = Rule(
    id="ABN_NUM_CHR_VALID",
    category="abnormality",
    description="Numerical abnormality chromosome must be 1-22, X, or Y",
    validate=_validate_numerical_chromosome,
)

breakpoint_arm_valid_rule = Rule(
    id="ABN_BP_ARM_VALID",
    category="abnormality",
    description="Breakpoint arm must be 'p' or 'q'",
    validate=_validate_breakpoint_arm,
)

inversion_two_breakpoints_rule = Rule(
    id="ABN_INV_TWO_BP",
    category="abnormality",
    description="Inversion must have exactly two breakpoints",
    validate=_validate_inversion_two_breakpoints,
)

translocation_breakpoint_count_rule = Rule(
    id="ABN_TRANS_BP_COUNT",
    category="abnormality",
    description="Translocation breakpoint count must match chromosome count",
    validate=_validate_translocation_breakpoint_count,
)

deletion_breakpoint_rule = Rule(
    id="ABN_DEL_BP",
    category="abnormality",
    description="Deletion must have 1-2 breakpoints, interstitial requires same arm",
    validate=_validate_deletion_breakpoints,
)

duplication_breakpoint_rule = Rule(
    id="ABN_DUP_BP",
    category="abnormality",
    description="Duplication must have 1-2 breakpoints, interstitial requires same arm",
    validate=_validate_duplication_breakpoints,
)

ring_chromosome_breakpoint_rule = Rule(
    id="ABN_RING_BP",
    category="abnormality",
    description="Ring chromosome must have 2 breakpoints on different arms",
    validate=_validate_ring_chromosome_breakpoints,
)

isochromosome_breakpoint_rule = Rule(
    id="ABN_ISO_BP",
    category="abnormality",
    description="Isochromosome must have exactly 1 breakpoint",
    validate=_validate_isochromosome_breakpoints,
)

triplication_breakpoint_rule = Rule(
    id="ABN_TRP_BP",
    category="abnormality",
    description="Triplication must have 2 breakpoints on same arm",
    validate=_validate_triplication_breakpoints,
)

quadruplication_breakpoint_rule = Rule(
    id="ABN_QDP_BP",
    category="abnormality",
    description="Quadruplication must have 2 breakpoints on same arm",
    validate=_validate_quadruplication_breakpoints,
)

dicentric_breakpoint_rule = Rule(
    id="ABN_DIC_BP",
    category="abnormality",
    description="Dicentric breakpoint count must match chromosome count",
    validate=_validate_dicentric_breakpoints,
)

isodicentric_breakpoint_rule = Rule(
    id="ABN_IDIC_BP",
    category="abnormality",
    description="Isodicentric must have exactly 1 breakpoint",
    validate=_validate_isodicentric_breakpoints,
)

robertsonian_breakpoint_rule = Rule(
    id="ABN_ROB_BP",
    category="abnormality",
    description="Robertsonian translocation breakpoint count must match chromosome count",
    validate=_validate_robertsonian_breakpoints,
)

add_breakpoint_rule = Rule(
    id="ABN_ADD_BP",
    category="abnormality",
    description="Add (additional material) must have exactly 1 breakpoint",
    validate=_validate_add_breakpoints,
)

fra_breakpoint_rule = Rule(
    id="ABN_FRA_BP",
    category="abnormality",
    description="Fragile site must have exactly 1 breakpoint",
    validate=_validate_fra_breakpoints,
)

ins_breakpoint_rule = Rule(
    id="ABN_INS_BP",
    category="abnormality",
    description="Insertion must have exactly 3 breakpoints",
    validate=_validate_ins_breakpoints,
)

dmin_breakpoint_rule = Rule(
    id="ABN_DMIN_BP",
    category="abnormality",
    description="Double minutes must have no breakpoints",
    validate=_validate_dmin_breakpoints,
)

hsr_breakpoint_rule = Rule(
    id="ABN_HSR_BP",
    category="abnormality",
    description="HSR must have 0 or 1 breakpoint",
    validate=_validate_hsr_breakpoints,
)

mar_breakpoint_rule = Rule(
    id="ABN_MAR_BP",
    category="abnormality",
    description="Marker chromosome must have no breakpoints",
    validate=_validate_mar_breakpoints,
)

pseudodicentric_breakpoint_rule = Rule(
    id="ABN_PSU_DIC_BP",
    category="abnormality",
    description="Pseudodicentric breakpoint count must match chromosome count",
    validate=_validate_pseudodicentric_breakpoints,
)

acentric_breakpoint_rule = Rule(
    id="ABN_ACE_BP",
    category="abnormality",
    description="Acentric fragment must have 1-2 breakpoints",
    validate=_validate_acentric_breakpoints,
)

telomeric_association_breakpoint_rule = Rule(
    id="ABN_TAS_BP",
    category="abnormality",
    description="Telomeric association breakpoint count must match chromosome count",
    validate=_validate_telomeric_association_breakpoints,
)

fission_breakpoint_rule = Rule(
    id="ABN_FIS_BP",
    category="abnormality",
    description="Fission must have exactly 1 breakpoint",
    validate=_validate_fission_breakpoints,
)

neocentromere_breakpoint_rule = Rule(
    id="ABN_NEO_BP",
    category="abnormality",
    description="Neocentromere must have exactly 1 breakpoint",
    validate=_validate_neocentromere_breakpoints,
)

incomplete_breakpoint_rule = Rule(
    id="ABN_INC_BP",
    category="abnormality",
    description="Incomplete karyotype marker must have no breakpoints",
    validate=_validate_incomplete_breakpoints,
)

# Export all rules
ALL_ABNORMALITY_RULES = [
    numerical_chromosome_valid_rule,
    breakpoint_arm_valid_rule,
    inversion_two_breakpoints_rule,
    translocation_breakpoint_count_rule,
    deletion_breakpoint_rule,
    duplication_breakpoint_rule,
    ring_chromosome_breakpoint_rule,
    isochromosome_breakpoint_rule,
    triplication_breakpoint_rule,
    quadruplication_breakpoint_rule,
    dicentric_breakpoint_rule,
    isodicentric_breakpoint_rule,
    robertsonian_breakpoint_rule,
    add_breakpoint_rule,
    fra_breakpoint_rule,
    ins_breakpoint_rule,
    dmin_breakpoint_rule,
    hsr_breakpoint_rule,
    mar_breakpoint_rule,
    pseudodicentric_breakpoint_rule,
    acentric_breakpoint_rule,
    telomeric_association_breakpoint_rule,
    fission_breakpoint_rule,
    neocentromere_breakpoint_rule,
    incomplete_breakpoint_rule,
]
