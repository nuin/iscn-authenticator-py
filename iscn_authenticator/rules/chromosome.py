# iscn_authenticator/rules/chromosome.py
"""Validation rules for chromosome count and sex chromosomes."""

from iscn_authenticator.models import KaryotypeAST
from iscn_authenticator.rules.base import Rule


def _validate_chromosome_count_numeric(ast: KaryotypeAST, _) -> list[str]:
    """Validate chromosome count is numeric."""
    if isinstance(ast.chromosome_count, str):
        # Range notation like "45~48" is valid
        if "~" in ast.chromosome_count:
            return []
        return [f"Chromosome count '{ast.chromosome_count}' is not numeric"]
    return []


def _validate_chromosome_count_range(ast: KaryotypeAST, _) -> list[str]:
    """Validate chromosome count is in valid range (23-92)."""
    count = ast.chromosome_count
    if isinstance(count, str):
        # Skip range validation for range notation
        return []
    if count < 23 or count > 92:
        return [f"Chromosome count {count} is outside valid range (must be between 23 and 92)"]
    return []


def _validate_sex_chromosomes_valid(ast: KaryotypeAST, _) -> list[str]:
    """Validate sex chromosomes contain at least one X."""
    sex = ast.sex_chromosomes
    if sex == "U":  # Undisclosed
        return []
    if "X" not in sex:
        return [f"Sex chromosomes '{sex}' must contain at least one X chromosome"]
    return []


def _validate_sex_chromosomes_coherence(ast: KaryotypeAST, _) -> list[str]:
    """Validate coherence between chromosome count and sex chromosome count."""
    count = ast.chromosome_count
    sex = ast.sex_chromosomes

    if isinstance(count, str) or sex == "U":
        # Skip coherence check for ranges or undisclosed
        return []

    sex_count = len(sex)

    # Basic coherence check for common karyotypes without listed abnormalities.
    # Only enforce strict coherence for 45 and 46 chromosome counts,
    # as these are the most common and well-defined cases.
    # Other counts (47+) often indicate aneuploidy that may not be fully
    # specified in the sex chromosomes alone.
    if not ast.abnormalities:
        if count == 46 and sex_count != 2:
            return [f"Chromosome count 46 requires 2 sex chromosomes, but found {sex_count} ('{sex}')"]
        if count == 45 and sex_count != 1:
            return [f"Chromosome count 45 requires 1 sex chromosome, but found {sex_count} ('{sex}')"]

    return []


# Rule instances
chromosome_count_numeric_rule = Rule(
    id="CHR_COUNT_NUMERIC",
    category="chromosome_count",
    description="Chromosome count must be numeric or valid range notation",
    validate=_validate_chromosome_count_numeric,
)

chromosome_count_range_rule = Rule(
    id="CHR_COUNT_RANGE",
    category="chromosome_count",
    description="Chromosome count must be between 23 and 92",
    validate=_validate_chromosome_count_range,
)

sex_chromosomes_valid_rule = Rule(
    id="SEX_CHR_VALID",
    category="sex_chromosomes",
    description="Sex chromosomes must contain at least one X",
    validate=_validate_sex_chromosomes_valid,
)

sex_chromosomes_coherence_rule = Rule(
    id="SEX_CHR_COHERENCE",
    category="coherence",
    description="Chromosome count must be coherent with sex chromosome count",
    validate=_validate_sex_chromosomes_coherence,
)

# Export all rules
ALL_CHROMOSOME_RULES = [
    chromosome_count_numeric_rule,
    chromosome_count_range_rule,
    sex_chromosomes_valid_rule,
    sex_chromosomes_coherence_rule,
]
