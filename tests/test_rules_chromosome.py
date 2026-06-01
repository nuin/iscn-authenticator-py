# tests/test_rules_chromosome.py
import unittest

from iscn_authenticator.models import KaryotypeAST
from iscn_authenticator.rules.chromosome import (
    chromosome_count_numeric_rule,
    chromosome_count_range_rule,
    sex_chromosomes_coherence_rule,
    sex_chromosomes_valid_rule,
)


class TestChromosomeCountRules(unittest.TestCase):
    def test_numeric_count_valid(self):
        ast = KaryotypeAST(46, "XX", [], None, None)
        errors = chromosome_count_numeric_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_numeric_count_valid_45(self):
        ast = KaryotypeAST(45, "X", [], None, None)
        errors = chromosome_count_numeric_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_numeric_count_range_notation_valid(self):
        ast = KaryotypeAST("45~48", "XX", [], None, None)
        errors = chromosome_count_numeric_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_count_range_valid(self):
        ast = KaryotypeAST(46, "XX", [], None, None)
        errors = chromosome_count_range_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_count_too_low(self):
        ast = KaryotypeAST(20, "XX", [], None, None)
        errors = chromosome_count_range_rule.validate(ast, ast)
        self.assertIn("between 23 and 92", errors[0])

    def test_count_too_high(self):
        ast = KaryotypeAST(100, "XX", [], None, None)
        errors = chromosome_count_range_rule.validate(ast, ast)
        self.assertIn("between 23 and 92", errors[0])


class TestSexChromosomeRules(unittest.TestCase):
    def test_valid_xx(self):
        ast = KaryotypeAST(46, "XX", [], None, None)
        errors = sex_chromosomes_valid_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_valid_xy(self):
        ast = KaryotypeAST(46, "XY", [], None, None)
        errors = sex_chromosomes_valid_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_valid_x(self):
        ast = KaryotypeAST(45, "X", [], None, None)
        errors = sex_chromosomes_valid_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_valid_xxy(self):
        ast = KaryotypeAST(47, "XXY", [], None, None)
        errors = sex_chromosomes_valid_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_valid_undisclosed(self):
        ast = KaryotypeAST(46, "U", [], None, None)
        errors = sex_chromosomes_valid_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_invalid_y_only(self):
        ast = KaryotypeAST(46, "Y", [], None, None)
        errors = sex_chromosomes_valid_rule.validate(ast, ast)
        self.assertIn("at least one X", errors[0])

    def test_invalid_yy(self):
        ast = KaryotypeAST(46, "YY", [], None, None)
        errors = sex_chromosomes_valid_rule.validate(ast, ast)
        self.assertIn("at least one X", errors[0])


class TestCoherenceRule(unittest.TestCase):
    def test_46_with_2_sex_chr(self):
        ast = KaryotypeAST(46, "XX", [], None, None)
        errors = sex_chromosomes_coherence_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_45_with_1_sex_chr(self):
        ast = KaryotypeAST(45, "X", [], None, None)
        errors = sex_chromosomes_coherence_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_47_with_3_sex_chr(self):
        ast = KaryotypeAST(47, "XXX", [], None, None)
        errors = sex_chromosomes_coherence_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_47_with_xxy(self):
        ast = KaryotypeAST(47, "XXY", [], None, None)
        errors = sex_chromosomes_coherence_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_47_with_xy_valid(self):
        # 47,XY could indicate +21 or other trisomy not listed in sex chr
        ast = KaryotypeAST(47, "XY", [], None, None)
        errors = sex_chromosomes_coherence_rule.validate(ast, ast)
        self.assertEqual(errors, [])

    def test_46_with_1_sex_chr_invalid(self):
        ast = KaryotypeAST(46, "X", [], None, None)
        errors = sex_chromosomes_coherence_rule.validate(ast, ast)
        self.assertIn("46", errors[0])
        self.assertIn("2", errors[0])

    def test_46_with_3_sex_chr_invalid(self):
        ast = KaryotypeAST(46, "XXY", [], None, None)
        errors = sex_chromosomes_coherence_rule.validate(ast, ast)
        self.assertIn("46", errors[0])

    def test_45_with_2_sex_chr_invalid(self):
        ast = KaryotypeAST(45, "XX", [], None, None)
        errors = sex_chromosomes_coherence_rule.validate(ast, ast)
        self.assertIn("45", errors[0])

    def test_range_notation_skips_coherence(self):
        ast = KaryotypeAST("45~48", "XX", [], None, None)
        errors = sex_chromosomes_coherence_rule.validate(ast, ast)
        self.assertEqual(errors, [])


if __name__ == "__main__":
    unittest.main()
