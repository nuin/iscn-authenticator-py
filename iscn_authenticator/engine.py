# iscn_authenticator/engine.py
"""Rule engine for applying validation rules to karyotype AST."""

from iscn_authenticator.models import KaryotypeAST, ValidationResult
from iscn_authenticator.rules.base import Rule


class RuleEngine:
    """Applies validation rules to a karyotype AST."""

    def __init__(self):
        self._rules: list[Rule] = []
        self._abnormality_rules: list[Rule] = []

    def add_rule(self, rule: Rule) -> None:
        """Add a rule to the engine."""
        self._rules.append(rule)

    def add_rules(self, rules: list[Rule]) -> None:
        """Add multiple rules to the engine."""
        self._rules.extend(rules)

    def add_abnormality_rule(self, rule: Rule) -> None:
        """Add an abnormality rule (applied to each abnormality)."""
        self._abnormality_rules.append(rule)

    def add_abnormality_rules(self, rules: list[Rule]) -> None:
        """Add multiple abnormality rules."""
        self._abnormality_rules.extend(rules)

    def validate(self, ast: KaryotypeAST) -> ValidationResult:
        """Validate a karyotype AST against all rules."""
        all_errors: list[str] = []

        # Apply AST-level rules
        for rule in self._rules:
            errors = rule.validate(ast, ast)
            all_errors.extend(errors)

        # Apply abnormality rules to each abnormality
        for abnormality in ast.abnormalities:
            for rule in self._abnormality_rules:
                errors = rule.validate(ast, abnormality)
                all_errors.extend(errors)

        return ValidationResult(valid=len(all_errors) == 0, errors=all_errors, parsed=ast)
