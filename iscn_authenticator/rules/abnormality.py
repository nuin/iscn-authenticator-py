# iscn_authenticator/rules/abnormality.py
"""Validation rules for karyotype abnormalities.

Validator implementations live in ``_abnormality_validators``; this module
exposes the public ``Rule`` instances and the ``ALL_ABNORMALITY_RULES`` list
used by the rule engine.
"""

from iscn_authenticator.rules import _abnormality_validators as v
from iscn_authenticator.rules.base import Rule

# Re-exported for backwards compatibility with prior versions of this module.
VALID_CHROMOSOMES = v.VALID_CHROMOSOMES

numerical_chromosome_valid_rule = Rule(
    id="ABN_NUM_CHR_VALID",
    category="abnormality",
    description="Numerical abnormality chromosome must be 1-22, X, or Y",
    validate=v.numerical_chromosome,
)

breakpoint_arm_valid_rule = Rule(
    id="ABN_BP_ARM_VALID",
    category="abnormality",
    description="Breakpoint arm must be 'p' or 'q'",
    validate=v.breakpoint_arm,
)

inversion_two_breakpoints_rule = Rule(
    id="ABN_INV_TWO_BP",
    category="abnormality",
    description="Inversion must have exactly two breakpoints",
    validate=v.inversion,
)

translocation_breakpoint_count_rule = Rule(
    id="ABN_TRANS_BP_COUNT",
    category="abnormality",
    description="Translocation breakpoint count must match chromosome count",
    validate=v.translocation,
)

deletion_breakpoint_rule = Rule(
    id="ABN_DEL_BP",
    category="abnormality",
    description="Deletion must have 1-2 breakpoints, interstitial requires same arm",
    validate=v.deletion,
)

duplication_breakpoint_rule = Rule(
    id="ABN_DUP_BP",
    category="abnormality",
    description="Duplication must have 1-2 breakpoints, interstitial requires same arm",
    validate=v.duplication,
)

ring_chromosome_breakpoint_rule = Rule(
    id="ABN_RING_BP",
    category="abnormality",
    description="Ring chromosome must have 2 breakpoints on different arms",
    validate=v.ring,
)

isochromosome_breakpoint_rule = Rule(
    id="ABN_ISO_BP",
    category="abnormality",
    description="Isochromosome must have exactly 1 breakpoint",
    validate=v.isochromosome,
)

triplication_breakpoint_rule = Rule(
    id="ABN_TRP_BP",
    category="abnormality",
    description="Triplication must have 2 breakpoints on same arm",
    validate=v.triplication,
)

quadruplication_breakpoint_rule = Rule(
    id="ABN_QDP_BP",
    category="abnormality",
    description="Quadruplication must have 2 breakpoints on same arm",
    validate=v.quadruplication,
)

dicentric_breakpoint_rule = Rule(
    id="ABN_DIC_BP",
    category="abnormality",
    description="Dicentric breakpoint count must match chromosome count",
    validate=v.dicentric,
)

isodicentric_breakpoint_rule = Rule(
    id="ABN_IDIC_BP",
    category="abnormality",
    description="Isodicentric must have exactly 1 breakpoint",
    validate=v.isodicentric,
)

robertsonian_breakpoint_rule = Rule(
    id="ABN_ROB_BP",
    category="abnormality",
    description="Robertsonian translocation breakpoint count must match chromosome count",
    validate=v.robertsonian,
)

add_breakpoint_rule = Rule(
    id="ABN_ADD_BP",
    category="abnormality",
    description="Add (additional material) must have exactly 1 breakpoint",
    validate=v.add_material,
)

fra_breakpoint_rule = Rule(
    id="ABN_FRA_BP",
    category="abnormality",
    description="Fragile site must have exactly 1 breakpoint",
    validate=v.fragile_site,
)

ins_breakpoint_rule = Rule(
    id="ABN_INS_BP",
    category="abnormality",
    description="Insertion must have exactly 3 breakpoints",
    validate=v.insertion,
)

dmin_breakpoint_rule = Rule(
    id="ABN_DMIN_BP",
    category="abnormality",
    description="Double minutes must have no breakpoints",
    validate=v.double_minutes,
)

hsr_breakpoint_rule = Rule(
    id="ABN_HSR_BP",
    category="abnormality",
    description="HSR must have 0 or 1 breakpoint",
    validate=v.hsr,
)

mar_breakpoint_rule = Rule(
    id="ABN_MAR_BP",
    category="abnormality",
    description="Marker chromosome must have no breakpoints",
    validate=v.marker,
)

pseudodicentric_breakpoint_rule = Rule(
    id="ABN_PSU_DIC_BP",
    category="abnormality",
    description="Pseudodicentric breakpoint count must match chromosome count",
    validate=v.pseudodicentric,
)

acentric_breakpoint_rule = Rule(
    id="ABN_ACE_BP",
    category="abnormality",
    description="Acentric fragment must have 1-2 breakpoints",
    validate=v.acentric,
)

telomeric_association_breakpoint_rule = Rule(
    id="ABN_TAS_BP",
    category="abnormality",
    description="Telomeric association breakpoint count must match chromosome count",
    validate=v.telomeric_assoc,
)

fission_breakpoint_rule = Rule(
    id="ABN_FIS_BP",
    category="abnormality",
    description="Fission must have exactly 1 breakpoint",
    validate=v.fission,
)

neocentromere_breakpoint_rule = Rule(
    id="ABN_NEO_BP",
    category="abnormality",
    description="Neocentromere must have exactly 1 breakpoint",
    validate=v.neocentromere,
)

incomplete_breakpoint_rule = Rule(
    id="ABN_INC_BP",
    category="abnormality",
    description="Incomplete karyotype marker must have no breakpoints",
    validate=v.incomplete,
)

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
