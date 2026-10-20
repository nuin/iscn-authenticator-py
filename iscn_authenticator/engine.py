# iscn_authenticator/engine.py
"""Rule engine for applying validation rules to karyotype AST."""

from iscn_authenticator._cell_line_context import cell_line_context
from iscn_authenticator.models import KaryotypeAST, ValidationResult
from iscn_authenticator.rules.base import Rule


class RuleEngine:
    """Applies validation rules to a karyotype AST."""

    def __init__(self):
        self._rules: list[Rule] = []
        self._abnormality_rules: list[Rule] = []
        self._karyotype_rules: list[Rule] = []

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
        errors = [error for rule in self._karyotype_rules for error in rule.validate(ast, ast)]
        if ast.cell_lines:
            for index, line in enumerate(ast.cell_lines, 1):
                context = cell_line_context(ast, line)
                errors.extend(f"Cell line {index}: {error}" for error in self._validate_line(context))
        else:
            errors.extend(self._validate_line(ast))
        return ValidationResult(valid=not errors, errors=errors, parsed=ast)

    def add_karyotype_rules(self, rules: list[Rule]) -> None:
        """Add rules evaluated once across the complete karyotype."""
        self._karyotype_rules.extend(rules)

    def _validate_line(self, ast: KaryotypeAST) -> list[str]:
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

        return all_errors
