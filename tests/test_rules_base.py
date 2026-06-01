# tests/test_rules_base.py
import unittest

from iscn_authenticator.rules.base import Rule


class TestRule(unittest.TestCase):
    def test_rule_creation(self):
        def validator(value, ast):
            return [] if value == 46 else ["Invalid count"]

        rule = Rule(
            id="CHR_COUNT_46",
            category="chromosome_count",
            description="Chromosome count must be 46",
            validate=validator,
        )
        self.assertEqual(rule.id, "CHR_COUNT_46")
        self.assertEqual(rule.category, "chromosome_count")

    def test_rule_validate_pass(self):
        def validator(value, ast):
            return [] if value == 46 else ["Invalid count"]

        rule = Rule(
            id="CHR_COUNT_46",
            category="chromosome_count",
            description="Test rule",
            validate=validator,
        )
        errors = rule.validate(46, None)
        self.assertEqual(errors, [])

    def test_rule_validate_fail(self):
        def validator(value, ast):
            return [] if value == 46 else [f"Expected 46, got {value}"]

        rule = Rule(
            id="CHR_COUNT_46",
            category="chromosome_count",
            description="Test rule",
            validate=validator,
        )
        errors = rule.validate(45, None)
        self.assertEqual(errors, ["Expected 46, got 45"])


if __name__ == "__main__":
    unittest.main()
