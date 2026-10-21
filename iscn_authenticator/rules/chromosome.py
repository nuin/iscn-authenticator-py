# iscn_authenticator/rules/chromosome.py
"""Validation rules for chromosome count and sex chromosomes."""

import re

from iscn_authenticator.models import KaryotypeAST
from iscn_authenticator.rules.base import Rule
from iscn_authenticator.rules.routine import numerical_total_rule


def _validate_chromosome_count_numeric(ast: KaryotypeAST, _) -> list[str]:
    """Validate chromosome count is numeric."""
    if isinstance(ast.chromosome_count, str):
        # Range notation like "45~48" is valid
        if "~" in ast.chromosome_count:
            match = re.fullmatch(r"([0-9]+)~([0-9]+)", ast.chromosome_count)
            if (
                match
                and max(len(match[1].lstrip("0")), len(match[2].lstrip("0"))) <= 16
                and 0 < int(match[1].lstrip("0") or "0") <= int(match[2].lstrip("0") or "0") <= 9_007_199_254_740_991
            ):
                return []
            return ["Chromosome count range must be ascending positive safe integers"]
        return [f"Chromosome count '{ast.chromosome_count}' is not numeric"]
    return []


def _validate_chromosome_count_range(ast: KaryotypeAST, _) -> list[str]:
    """Validate chromosome count is is a positive safe integer."""
    count = ast.chromosome_count
    if count is None or isinstance(count, str):
        # Skip range validation for range notation
        return []
    if type(count) is not int or not 0 < count <= 9_007_199_254_740_991:
        return [f"Chromosome count {count} must be a positive safe integer"]
    return []


def _validate_sex_chromosomes_valid(ast: KaryotypeAST, _) -> list[str]:
    """Validate sex chromosomes contain at least one X."""
    sex = ast.sex_chromosomes
    if ast.chromosome_count is None or sex == "U":  # Undisclosed
        return []
    if sex == "Y" and any(abn.type == "-" and abn.chromosome == "X" for abn in ast.abnormalities):
        return []
    if "X" not in sex:
        return [f"Sex chromosomes '{sex}' must contain at least one X chromosome"]
    return []


def _validate_sex_chromosomes_coherence(ast: KaryotypeAST, _) -> list[str]:
    """Validate coherence between chromosome count and sex chromosome count."""
    count = ast.chromosome_count
    sex = ast.sex_chromosomes

    if count is None or isinstance(count, str) or sex == "U":
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
    description="Chromosome count must be a positive safe integer",
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
    numerical_total_rule,
]
