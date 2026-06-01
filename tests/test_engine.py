# tests/test_engine.py
import unittest

from iscn_authenticator.engine import RuleEngine
from iscn_authenticator.models import KaryotypeAST
from iscn_authenticator.rules.base import Rule


class TestRuleEngine(unittest.TestCase):
    def test_engine_no_rules(self):
        engine = RuleEngine()
        ast = KaryotypeAST(
            chromosome_count=46,
            sex_chromosomes="XX",
            abnormalities=[],
            cell_lines=None,
            modifiers=None,
        )
        result = engine.validate(ast)
        self.assertTrue(result.valid)
        self.assertEqual(result.errors, [])

    def test_engine_with_passing_rule(self):
        engine = RuleEngine()
        engine.add_rule(
            Rule(
                id="TEST_PASS",
                category="test",
                description="Always passes",
                validate=lambda v, ast: [],
            )
        )
        ast = KaryotypeAST(
            chromosome_count=46,
            sex_chromosomes="XX",
            abnormalities=[],
            cell_lines=None,
            modifiers=None,
        )
        result = engine.validate(ast)
        self.assertTrue(result.valid)

    def test_engine_with_failing_rule(self):
        engine = RuleEngine()
        engine.add_rule(
            Rule(
                id="TEST_FAIL",
                category="test",
                description="Always fails",
                validate=lambda v, ast: ["Test error"],
            )
        )
        ast = KaryotypeAST(
            chromosome_count=46,
            sex_chromosomes="XX",
            abnormalities=[],
            cell_lines=None,
            modifiers=None,
        )
        result = engine.validate(ast)
        self.assertFalse(result.valid)
        self.assertIn("Test error", result.errors)

    def test_engine_collects_all_errors(self):
        engine = RuleEngine()
        engine.add_rule(
            Rule(
                id="FAIL_1",
                category="test",
                description="Fails with error 1",
                validate=lambda v, ast: ["Error 1"],
            )
        )
        engine.add_rule(
            Rule(
                id="FAIL_2",
                category="test",
                description="Fails with error 2",
                validate=lambda v, ast: ["Error 2"],
            )
        )
        ast = KaryotypeAST(
            chromosome_count=46,
            sex_chromosomes="XX",
            abnormalities=[],
            cell_lines=None,
            modifiers=None,
        )
        result = engine.validate(ast)
        self.assertFalse(result.valid)
        self.assertEqual(len(result.errors), 2)
        self.assertIn("Error 1", result.errors)
        self.assertIn("Error 2", result.errors)

    def test_engine_add_rules_multiple(self):
        engine = RuleEngine()
        rules = [
            Rule(
                id="RULE_1",
                category="test",
                description="Rule 1",
                validate=lambda v, ast: [],
            ),
            Rule(
                id="RULE_2",
                category="test",
                description="Rule 2",
                validate=lambda v, ast: [],
            ),
        ]
        engine.add_rules(rules)
        ast = KaryotypeAST(46, "XX", [], None, None)
        result = engine.validate(ast)
        self.assertTrue(result.valid)

    def test_engine_result_contains_parsed_ast(self):
        engine = RuleEngine()
        ast = KaryotypeAST(
            chromosome_count=47,
            sex_chromosomes="XXY",
            abnormalities=[],
            cell_lines=None,
            modifiers=None,
        )
        result = engine.validate(ast)
        self.assertEqual(result.parsed, ast)


if __name__ == "__main__":
    unittest.main()
