# iscn_authenticator/main.py
"""ISCN Karyotype Validation API."""

import re

from iscn_authenticator.engine import RuleEngine
from iscn_authenticator.explain import explain
from iscn_authenticator.models import ValidationResult
from iscn_authenticator.parser import KaryotypeParser, ParseError
from iscn_authenticator.rules.abnormality import ALL_ABNORMALITY_RULES
from iscn_authenticator.rules.chromosome import ALL_CHROMOSOME_RULES

_parser = KaryotypeParser()
_engine = RuleEngine()
_engine.add_rules(ALL_CHROMOSOME_RULES)
_engine.add_abnormality_rules(ALL_ABNORMALITY_RULES)


def validate_karyotype(karyotype: str) -> ValidationResult:
    """
    Validate an ISCN karyotype string.

    Returns a ValidationResult containing:
    - valid: bool indicating if karyotype is valid
    - errors: list of error messages if invalid
    - parsed: KaryotypeAST if parsing succeeded
    - explanation: Optional[ExplainResult] providing clinical context
    """
    try:
        ast = _parser.parse(karyotype)
    except ParseError as e:
        return ValidationResult(valid=False, errors=[str(e)], parsed=None)

    result = _engine.validate(ast)

    # Populate explanations
    if result.parsed:
        result.explanation = explain(result.parsed)
        for abn in result.parsed.abnormalities:
            abn.explanation = explain(abn)

    return result


def is_valid_karyotype(karyotype: str) -> bool:
    """
    Validates a simple ISCN karyotype string.

    Backward-compatible API that returns only True/False.
    For detailed errors, use validate_karyotype() instead.
    """
    return validate_karyotype(karyotype).valid


# Legacy functions kept for compatibility with existing code
# These can be removed in a future version


def _validate_total_chromosome_number(number_part: str) -> bool:
    """Validates the total chromosome number part of the karyotype."""
    return number_part.isdigit()


def _validate_sex_chromosomes(sex_chromosome_part: str) -> bool:
    """Validates the sex chromosomes part of the karyotype."""
    return bool(re.match(r"^(X|Y)+$", sex_chromosome_part))


def _validate_coherence(total_chromosome_number: int, sex_chromosomes: str) -> bool:
    """Validates coherence between total chromosome number and sex chromosomes."""
    if total_chromosome_number == 46:
        if len(sex_chromosomes) != 2:
            return False
    elif total_chromosome_number == 45 and len(sex_chromosomes) != 1:
        return False
    return True


def _validate_deletion_content(content: str) -> bool:
    """Validates the content string of a deletion, e.g., 'q13' or 'q13q33'."""
    # Terminal deletion: [pq]\d+(\.\d+)?
    if re.fullmatch(r"[pq]\d+(\.\d+)?", content):
        return True

    # Interstitial deletion: [pq]\d+(\.\d+)?[pq]\d+(\.\d+)?
    return bool(re.fullmatch("[pq]\\d+(\\.\\d+)?[pq]\\d+(\\.\\d+)?", content))


def _validate_abnormalities(abnormality_parts: list[str]) -> bool:
    """Validates the abnormality parts of the karyotype."""
    for part in abnormality_parts:
        if re.fullmatch(r"^[+-]\d+$", part):
            continue  # It is a valid numeric abnormality

        del_match = re.fullmatch(r"^del\((\d{1,2}|[XY])\)\((.+)\)$", part)
        if del_match:
            content = del_match.group(2)
            if _validate_deletion_content(content):
                continue

        return False  # If no match, then it's an invalid part
    return True


if __name__ == "__main__":
    karyotype_string = input("Enter karyotype string: ")
    result = validate_karyotype(karyotype_string)
    if result.valid:
        print("Valid")
    else:
        print("Invalid")
        for error in result.errors:
            print(f"  - {error}")
