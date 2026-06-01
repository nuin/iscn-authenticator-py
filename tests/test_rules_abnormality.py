# tests/test_rules_abnormality.py
import unittest

from iscn_authenticator.models import Abnormality, Breakpoint, KaryotypeAST
from iscn_authenticator.rules.abnormality import (
    breakpoint_arm_valid_rule,
    inversion_two_breakpoints_rule,
    numerical_chromosome_valid_rule,
    translocation_breakpoint_count_rule,
)


class TestNumericalChromosomeRule(unittest.TestCase):
    def test_valid_autosome_gain(self):
        abn = Abnormality("+", "21", [], None, False, None, "+21")
        ast = KaryotypeAST(47, "XX", [abn], None, None)
        errors = numerical_chromosome_valid_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_valid_autosome_loss(self):
        abn = Abnormality("-", "7", [], None, False, None, "-7")
        ast = KaryotypeAST(45, "XY", [abn], None, None)
        errors = numerical_chromosome_valid_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_valid_sex_chromosome_x(self):
        abn = Abnormality("+", "X", [], None, False, None, "+X")
        ast = KaryotypeAST(47, "XXX", [abn], None, None)
        errors = numerical_chromosome_valid_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_valid_sex_chromosome_y(self):
        abn = Abnormality("-", "Y", [], None, False, None, "-Y")
        ast = KaryotypeAST(45, "X", [abn], None, None)
        errors = numerical_chromosome_valid_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_chromosome_0(self):
        abn = Abnormality("+", "0", [], None, False, None, "+0")
        ast = KaryotypeAST(47, "XX", [abn], None, None)
        errors = numerical_chromosome_valid_rule.validate(ast, abn)
        self.assertIn("Invalid chromosome", errors[0])

    def test_invalid_chromosome_23(self):
        abn = Abnormality("+", "23", [], None, False, None, "+23")
        ast = KaryotypeAST(47, "XX", [abn], None, None)
        errors = numerical_chromosome_valid_rule.validate(ast, abn)
        self.assertIn("Invalid chromosome", errors[0])

    def test_invalid_chromosome_99(self):
        abn = Abnormality("+", "99", [], None, False, None, "+99")
        ast = KaryotypeAST(47, "XX", [abn], None, None)
        errors = numerical_chromosome_valid_rule.validate(ast, abn)
        self.assertIn("Invalid chromosome", errors[0])

    def test_skips_non_numerical(self):
        abn = Abnormality("del", "5", [], None, False, None, "del(5)(q13)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = numerical_chromosome_valid_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestBreakpointArmRule(unittest.TestCase):
    def test_valid_p_arm(self):
        bp = Breakpoint("p", 1, 3, None, False)
        abn = Abnormality("del", "5", [bp], None, False, None, "del(5)(p13)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = breakpoint_arm_valid_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_valid_q_arm(self):
        bp = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("del", "5", [bp], None, False, None, "del(5)(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = breakpoint_arm_valid_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_skips_numerical_abnormality(self):
        abn = Abnormality("+", "21", [], None, False, None, "+21")
        ast = KaryotypeAST(47, "XX", [abn], None, None)
        errors = breakpoint_arm_valid_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestInversionTwoBreakpointsRule(unittest.TestCase):
    def test_valid_two_breakpoints(self):
        bp1 = Breakpoint("p", 1, 2, None, False)
        bp2 = Breakpoint("q", 1, 3, None, False)
        abn = Abnormality("inv", "9", [bp1, bp2], None, False, None, "inv(9)(p12q13)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = inversion_two_breakpoints_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_one_breakpoint(self):
        bp1 = Breakpoint("p", 1, 2, None, False)
        abn = Abnormality("inv", "9", [bp1], None, False, None, "inv(9)(p12)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = inversion_two_breakpoints_rule.validate(ast, abn)
        self.assertIn("two breakpoints", errors[0].lower())

    def test_invalid_three_breakpoints(self):
        bp1 = Breakpoint("p", 1, 2, None, False)
        bp2 = Breakpoint("q", 1, 3, None, False)
        bp3 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("inv", "9", [bp1, bp2, bp3], None, False, None, "inv(9)(p12q13q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = inversion_two_breakpoints_rule.validate(ast, abn)
        self.assertIn("two breakpoints", errors[0].lower())

    def test_skips_non_inversion(self):
        bp1 = Breakpoint("q", 1, 3, None, False)
        abn = Abnormality("del", "5", [bp1], None, False, None, "del(5)(q13)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = inversion_two_breakpoints_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestDeletionBreakpointRule(unittest.TestCase):
    def test_valid_terminal_deletion_one_breakpoint(self):
        """Terminal deletion has one breakpoint."""
        from iscn_authenticator.rules.abnormality import deletion_breakpoint_rule

        bp1 = Breakpoint("q", 1, 3, None, False)
        abn = Abnormality("del", "5", [bp1], None, False, None, "del(5)(q13)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = deletion_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_valid_interstitial_deletion_two_breakpoints_same_arm(self):
        """Interstitial deletion has two breakpoints on same arm."""
        from iscn_authenticator.rules.abnormality import deletion_breakpoint_rule

        bp1 = Breakpoint("q", 1, 3, None, False)
        bp2 = Breakpoint("q", 3, 3, None, False)
        abn = Abnormality("del", "5", [bp1, bp2], None, False, None, "del(5)(q13q33)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = deletion_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_interstitial_deletion_different_arms(self):
        """Interstitial deletion cannot have breakpoints on different arms."""
        from iscn_authenticator.rules.abnormality import deletion_breakpoint_rule

        bp1 = Breakpoint("p", 1, 2, None, False)
        bp2 = Breakpoint("q", 1, 3, None, False)
        abn = Abnormality("del", "5", [bp1, bp2], None, False, None, "del(5)(p12q13)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = deletion_breakpoint_rule.validate(ast, abn)
        self.assertIn("same arm", errors[0].lower())

    def test_invalid_deletion_three_breakpoints(self):
        """Deletion cannot have three breakpoints."""
        from iscn_authenticator.rules.abnormality import deletion_breakpoint_rule

        bp1 = Breakpoint("q", 1, 3, None, False)
        bp2 = Breakpoint("q", 2, 1, None, False)
        bp3 = Breakpoint("q", 3, 3, None, False)
        abn = Abnormality("del", "5", [bp1, bp2, bp3], None, False, None, "del(5)(q13q21q33)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = deletion_breakpoint_rule.validate(ast, abn)
        self.assertIn("one or two breakpoints", errors[0].lower())

    def test_skips_non_deletion(self):
        """Rule only applies to deletions."""
        from iscn_authenticator.rules.abnormality import deletion_breakpoint_rule

        bp1 = Breakpoint("q", 1, 3, None, False)
        bp2 = Breakpoint("q", 3, 3, None, False)
        abn = Abnormality("dup", "5", [bp1, bp2], None, False, None, "dup(5)(q13q33)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = deletion_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestTranslocationBreakpointCountRule(unittest.TestCase):
    def test_valid_two_chromosome_two_breakpoint(self):
        bp1 = Breakpoint("q", 3, 4, None, False)
        bp2 = Breakpoint("q", 1, 1, "2", False)
        abn = Abnormality("t", "9;22", [bp1, bp2], None, False, None, "t(9;22)(q34;q11.2)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = translocation_breakpoint_count_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_valid_three_chromosome_three_breakpoint(self):
        bp1 = Breakpoint("p", 3, 2, None, False)
        bp2 = Breakpoint("q", 2, 1, None, False)
        bp3 = Breakpoint("q", 3, 1, None, False)
        abn = Abnormality("t", "1;3;5", [bp1, bp2, bp3], None, False, None, "t(1;3;5)(p32;q21;q31)")
        ast = KaryotypeAST(46, "XY", [abn], None, None)
        errors = translocation_breakpoint_count_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_two_chromosomes_one_breakpoint(self):
        bp1 = Breakpoint("q", 3, 4, None, False)
        abn = Abnormality("t", "9;22", [bp1], None, False, None, "t(9;22)(q34)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = translocation_breakpoint_count_rule.validate(ast, abn)
        self.assertIn("breakpoints", errors[0].lower())

    def test_invalid_three_chromosomes_two_breakpoints(self):
        bp1 = Breakpoint("p", 3, 2, None, False)
        bp2 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("t", "1;3;5", [bp1, bp2], None, False, None, "t(1;3;5)(p32;q21)")
        ast = KaryotypeAST(46, "XY", [abn], None, None)
        errors = translocation_breakpoint_count_rule.validate(ast, abn)
        self.assertIn("breakpoints", errors[0].lower())

    def test_skips_non_translocation(self):
        bp1 = Breakpoint("p", 1, 2, None, False)
        bp2 = Breakpoint("q", 1, 3, None, False)
        abn = Abnormality("inv", "9", [bp1, bp2], None, False, None, "inv(9)(p12q13)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = translocation_breakpoint_count_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestAddBreakpointRule(unittest.TestCase):
    def test_valid_add_one_breakpoint(self):
        """Add with one breakpoint."""
        from iscn_authenticator.rules.abnormality import add_breakpoint_rule

        bp1 = Breakpoint("p", 2, 2, None, False)
        abn = Abnormality("add", "7", [bp1], None, False, None, "add(7)(p22)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = add_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_add_no_breakpoint(self):
        """Add must have a breakpoint."""
        from iscn_authenticator.rules.abnormality import add_breakpoint_rule

        abn = Abnormality("add", "7", [], None, False, None, "add(7)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = add_breakpoint_rule.validate(ast, abn)
        self.assertIn("one breakpoint", errors[0].lower())

    def test_invalid_add_two_breakpoints(self):
        """Add cannot have two breakpoints."""
        from iscn_authenticator.rules.abnormality import add_breakpoint_rule

        bp1 = Breakpoint("p", 2, 2, None, False)
        bp2 = Breakpoint("q", 1, 1, None, False)
        abn = Abnormality("add", "7", [bp1, bp2], None, False, None, "add(7)(p22q11)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = add_breakpoint_rule.validate(ast, abn)
        self.assertIn("one breakpoint", errors[0].lower())

    def test_skips_non_add(self):
        """Rule only applies to add."""
        from iscn_authenticator.rules.abnormality import add_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("del", "5", [bp1], None, False, None, "del(5)(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = add_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestRobertsonianBreakpointRule(unittest.TestCase):
    def test_valid_robertsonian_two_chromosomes_two_breakpoints(self):
        """Robertsonian with two chromosomes and two breakpoints."""
        from iscn_authenticator.rules.abnormality import robertsonian_breakpoint_rule

        bp1 = Breakpoint("q", 1, 0, None, False)
        bp2 = Breakpoint("q", 1, 0, None, False)
        abn = Abnormality("rob", "13;14", [bp1, bp2], None, False, None, "rob(13;14)(q10;q10)")
        ast = KaryotypeAST(45, "XX", [abn], None, None)
        errors = robertsonian_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_robertsonian_two_chromosomes_one_breakpoint(self):
        """Robertsonian with two chromosomes but one breakpoint."""
        from iscn_authenticator.rules.abnormality import robertsonian_breakpoint_rule

        bp1 = Breakpoint("q", 1, 0, None, False)
        abn = Abnormality("rob", "13;14", [bp1], None, False, None, "rob(13;14)(q10)")
        ast = KaryotypeAST(45, "XX", [abn], None, None)
        errors = robertsonian_breakpoint_rule.validate(ast, abn)
        self.assertIn("breakpoints", errors[0].lower())

    def test_invalid_robertsonian_two_chromosomes_three_breakpoints(self):
        """Robertsonian with two chromosomes but three breakpoints."""
        from iscn_authenticator.rules.abnormality import robertsonian_breakpoint_rule

        bp1 = Breakpoint("q", 1, 0, None, False)
        bp2 = Breakpoint("q", 1, 0, None, False)
        bp3 = Breakpoint("p", 1, 0, None, False)
        abn = Abnormality(
            "rob",
            "13;14",
            [bp1, bp2, bp3],
            None,
            False,
            None,
            "rob(13;14)(q10;q10;p10)",
        )
        ast = KaryotypeAST(45, "XX", [abn], None, None)
        errors = robertsonian_breakpoint_rule.validate(ast, abn)
        self.assertIn("breakpoints", errors[0].lower())

    def test_skips_non_robertsonian(self):
        """Rule only applies to Robertsonian translocations."""
        from iscn_authenticator.rules.abnormality import robertsonian_breakpoint_rule

        bp1 = Breakpoint("q", 3, 4, None, False)
        bp2 = Breakpoint("q", 1, 1, None, False)
        abn = Abnormality("t", "9;22", [bp1, bp2], None, False, None, "t(9;22)(q34;q11)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = robertsonian_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestIsodicentricBreakpointRule(unittest.TestCase):
    def test_valid_isodicentric_one_breakpoint(self):
        """Isodicentric with one breakpoint."""
        from iscn_authenticator.rules.abnormality import isodicentric_breakpoint_rule

        bp1 = Breakpoint("q", 1, 1, None, False)
        abn = Abnormality("idic", "Y", [bp1], None, False, None, "idic(Y)(q11)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = isodicentric_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_isodicentric_no_breakpoint(self):
        """Isodicentric must have a breakpoint."""
        from iscn_authenticator.rules.abnormality import isodicentric_breakpoint_rule

        abn = Abnormality("idic", "Y", [], None, False, None, "idic(Y)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = isodicentric_breakpoint_rule.validate(ast, abn)
        self.assertIn("one breakpoint", errors[0].lower())

    def test_invalid_isodicentric_two_breakpoints(self):
        """Isodicentric cannot have two breakpoints."""
        from iscn_authenticator.rules.abnormality import isodicentric_breakpoint_rule

        bp1 = Breakpoint("q", 1, 1, None, False)
        bp2 = Breakpoint("p", 1, 1, None, False)
        abn = Abnormality("idic", "Y", [bp1, bp2], None, False, None, "idic(Y)(q11p11)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = isodicentric_breakpoint_rule.validate(ast, abn)
        self.assertIn("one breakpoint", errors[0].lower())

    def test_skips_non_isodicentric(self):
        """Rule only applies to isodicentric chromosomes."""
        from iscn_authenticator.rules.abnormality import isodicentric_breakpoint_rule

        bp1 = Breakpoint("q", 1, 0, None, False)
        abn = Abnormality("i", "17", [bp1], None, False, None, "i(17)(q10)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = isodicentric_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestDicentricBreakpointRule(unittest.TestCase):
    def test_valid_dicentric_two_chromosomes_two_breakpoints(self):
        """Dicentric with two chromosomes and two breakpoints."""
        from iscn_authenticator.rules.abnormality import dicentric_breakpoint_rule

        bp1 = Breakpoint("q", 1, 4, None, False)
        bp2 = Breakpoint("q", 1, 1, None, False)
        abn = Abnormality("dic", "13;14", [bp1, bp2], None, False, None, "dic(13;14)(q14;q11)")
        ast = KaryotypeAST(45, "XX", [abn], None, None)
        errors = dicentric_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_dicentric_two_chromosomes_one_breakpoint(self):
        """Dicentric with two chromosomes but one breakpoint."""
        from iscn_authenticator.rules.abnormality import dicentric_breakpoint_rule

        bp1 = Breakpoint("q", 1, 4, None, False)
        abn = Abnormality("dic", "13;14", [bp1], None, False, None, "dic(13;14)(q14)")
        ast = KaryotypeAST(45, "XX", [abn], None, None)
        errors = dicentric_breakpoint_rule.validate(ast, abn)
        self.assertIn("breakpoints", errors[0].lower())

    def test_invalid_dicentric_two_chromosomes_three_breakpoints(self):
        """Dicentric with two chromosomes but three breakpoints."""
        from iscn_authenticator.rules.abnormality import dicentric_breakpoint_rule

        bp1 = Breakpoint("q", 1, 4, None, False)
        bp2 = Breakpoint("q", 1, 1, None, False)
        bp3 = Breakpoint("p", 1, 2, None, False)
        abn = Abnormality(
            "dic",
            "13;14",
            [bp1, bp2, bp3],
            None,
            False,
            None,
            "dic(13;14)(q14;q11;p12)",
        )
        ast = KaryotypeAST(45, "XX", [abn], None, None)
        errors = dicentric_breakpoint_rule.validate(ast, abn)
        self.assertIn("breakpoints", errors[0].lower())

    def test_skips_non_dicentric(self):
        """Rule only applies to dicentric chromosomes."""
        from iscn_authenticator.rules.abnormality import dicentric_breakpoint_rule

        bp1 = Breakpoint("q", 3, 4, None, False)
        bp2 = Breakpoint("q", 1, 1, None, False)
        abn = Abnormality("t", "9;22", [bp1, bp2], None, False, None, "t(9;22)(q34;q11)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = dicentric_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestQuadruplicationBreakpointRule(unittest.TestCase):
    def test_valid_quadruplication_two_breakpoints_same_arm(self):
        """Quadruplication with two breakpoints on same arm."""
        from iscn_authenticator.rules.abnormality import quadruplication_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 3, 2, None, False)
        abn = Abnormality("qdp", "1", [bp1, bp2], None, False, None, "qdp(1)(q21q32)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = quadruplication_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_quadruplication_one_breakpoint(self):
        """Quadruplication should have two breakpoints."""
        from iscn_authenticator.rules.abnormality import quadruplication_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("qdp", "1", [bp1], None, False, None, "qdp(1)(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = quadruplication_breakpoint_rule.validate(ast, abn)
        self.assertIn("two breakpoints", errors[0].lower())

    def test_invalid_quadruplication_different_arms(self):
        """Quadruplication breakpoints must be on same arm."""
        from iscn_authenticator.rules.abnormality import quadruplication_breakpoint_rule

        bp1 = Breakpoint("p", 1, 3, None, False)
        bp2 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("qdp", "1", [bp1, bp2], None, False, None, "qdp(1)(p13q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = quadruplication_breakpoint_rule.validate(ast, abn)
        self.assertIn("same arm", errors[0].lower())

    def test_skips_non_quadruplication(self):
        """Rule only applies to quadruplications."""
        from iscn_authenticator.rules.abnormality import quadruplication_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 3, 1, None, False)
        abn = Abnormality("trp", "1", [bp1, bp2], None, False, None, "trp(1)(q21q31)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = quadruplication_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestTriplicationBreakpointRule(unittest.TestCase):
    def test_valid_triplication_two_breakpoints_same_arm(self):
        """Triplication with two breakpoints on same arm."""
        from iscn_authenticator.rules.abnormality import triplication_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 3, 2, None, False)
        abn = Abnormality("trp", "1", [bp1, bp2], None, False, None, "trp(1)(q21q32)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = triplication_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_triplication_one_breakpoint(self):
        """Triplication should have two breakpoints."""
        from iscn_authenticator.rules.abnormality import triplication_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("trp", "1", [bp1], None, False, None, "trp(1)(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = triplication_breakpoint_rule.validate(ast, abn)
        self.assertIn("two breakpoints", errors[0].lower())

    def test_invalid_triplication_different_arms(self):
        """Triplication breakpoints must be on same arm."""
        from iscn_authenticator.rules.abnormality import triplication_breakpoint_rule

        bp1 = Breakpoint("p", 1, 3, None, False)
        bp2 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("trp", "1", [bp1, bp2], None, False, None, "trp(1)(p13q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = triplication_breakpoint_rule.validate(ast, abn)
        self.assertIn("same arm", errors[0].lower())

    def test_invalid_triplication_three_breakpoints(self):
        """Triplication cannot have three breakpoints."""
        from iscn_authenticator.rules.abnormality import triplication_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 2, 5, None, False)
        bp3 = Breakpoint("q", 3, 2, None, False)
        abn = Abnormality("trp", "1", [bp1, bp2, bp3], None, False, None, "trp(1)(q21q25q32)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = triplication_breakpoint_rule.validate(ast, abn)
        self.assertIn("two breakpoints", errors[0].lower())

    def test_skips_non_triplication(self):
        """Rule only applies to triplications."""
        from iscn_authenticator.rules.abnormality import triplication_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 3, 1, None, False)
        abn = Abnormality("dup", "1", [bp1, bp2], None, False, None, "dup(1)(q21q31)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = triplication_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestIsochromosomeBreakpointRule(unittest.TestCase):
    def test_valid_isochromosome_q10(self):
        """Isochromosome with q10 breakpoint (centromere)."""
        from iscn_authenticator.rules.abnormality import isochromosome_breakpoint_rule

        bp1 = Breakpoint("q", 1, 0, None, False)
        abn = Abnormality("i", "17", [bp1], None, False, None, "i(17)(q10)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = isochromosome_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_valid_isochromosome_p10(self):
        """Isochromosome with p10 breakpoint (centromere)."""
        from iscn_authenticator.rules.abnormality import isochromosome_breakpoint_rule

        bp1 = Breakpoint("p", 1, 0, None, False)
        abn = Abnormality("i", "17", [bp1], None, False, None, "i(17)(p10)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = isochromosome_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_isochromosome_no_breakpoint(self):
        """Isochromosome must have a breakpoint."""
        from iscn_authenticator.rules.abnormality import isochromosome_breakpoint_rule

        abn = Abnormality("i", "17", [], None, False, None, "i(17)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = isochromosome_breakpoint_rule.validate(ast, abn)
        self.assertIn("one breakpoint", errors[0].lower())

    def test_invalid_isochromosome_two_breakpoints(self):
        """Isochromosome cannot have two breakpoints."""
        from iscn_authenticator.rules.abnormality import isochromosome_breakpoint_rule

        bp1 = Breakpoint("q", 1, 0, None, False)
        bp2 = Breakpoint("p", 1, 0, None, False)
        abn = Abnormality("i", "17", [bp1, bp2], None, False, None, "i(17)(q10p10)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = isochromosome_breakpoint_rule.validate(ast, abn)
        self.assertIn("one breakpoint", errors[0].lower())

    def test_skips_non_isochromosome(self):
        """Rule only applies to isochromosomes."""
        from iscn_authenticator.rules.abnormality import isochromosome_breakpoint_rule

        bp1 = Breakpoint("p", 1, 2, None, False)
        bp2 = Breakpoint("q", 1, 3, None, False)
        abn = Abnormality("inv", "9", [bp1, bp2], None, False, None, "inv(9)(p12q13)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = isochromosome_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestRingChromosomeBreakpointRule(unittest.TestCase):
    def test_valid_ring_two_breakpoints_different_arms(self):
        """Ring chromosome requires two breakpoints on different arms."""
        from iscn_authenticator.rules.abnormality import ring_chromosome_breakpoint_rule

        bp1 = Breakpoint("p", 1, 1, None, False)
        bp2 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("r", "7", [bp1, bp2], None, False, None, "r(7)(p11q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = ring_chromosome_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_ring_one_breakpoint(self):
        """Ring chromosome cannot have one breakpoint."""
        from iscn_authenticator.rules.abnormality import ring_chromosome_breakpoint_rule

        bp1 = Breakpoint("p", 1, 1, None, False)
        abn = Abnormality("r", "7", [bp1], None, False, None, "r(7)(p11)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = ring_chromosome_breakpoint_rule.validate(ast, abn)
        self.assertIn("two breakpoints", errors[0].lower())

    def test_invalid_ring_same_arm(self):
        """Ring chromosome breakpoints must be on different arms."""
        from iscn_authenticator.rules.abnormality import ring_chromosome_breakpoint_rule

        bp1 = Breakpoint("p", 1, 1, None, False)
        bp2 = Breakpoint("p", 1, 3, None, False)
        abn = Abnormality("r", "7", [bp1, bp2], None, False, None, "r(7)(p11p13)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = ring_chromosome_breakpoint_rule.validate(ast, abn)
        self.assertIn("different arms", errors[0].lower())

    def test_invalid_ring_three_breakpoints(self):
        """Ring chromosome cannot have three breakpoints."""
        from iscn_authenticator.rules.abnormality import ring_chromosome_breakpoint_rule

        bp1 = Breakpoint("p", 1, 1, None, False)
        bp2 = Breakpoint("q", 2, 1, None, False)
        bp3 = Breakpoint("q", 3, 1, None, False)
        abn = Abnormality("r", "7", [bp1, bp2, bp3], None, False, None, "r(7)(p11q21q31)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = ring_chromosome_breakpoint_rule.validate(ast, abn)
        self.assertIn("two breakpoints", errors[0].lower())

    def test_skips_non_ring(self):
        """Rule only applies to ring chromosomes."""
        from iscn_authenticator.rules.abnormality import ring_chromosome_breakpoint_rule

        bp1 = Breakpoint("p", 1, 2, None, False)
        bp2 = Breakpoint("q", 1, 3, None, False)
        abn = Abnormality("inv", "9", [bp1, bp2], None, False, None, "inv(9)(p12q13)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = ring_chromosome_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestDuplicationBreakpointRule(unittest.TestCase):
    def test_valid_tandem_duplication_one_breakpoint(self):
        """Tandem duplication has one breakpoint."""
        from iscn_authenticator.rules.abnormality import duplication_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("dup", "1", [bp1], None, False, None, "dup(1)(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = duplication_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_valid_duplication_two_breakpoints_same_arm(self):
        """Duplication with two breakpoints on same arm."""
        from iscn_authenticator.rules.abnormality import duplication_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 3, 1, None, False)
        abn = Abnormality("dup", "1", [bp1, bp2], None, False, None, "dup(1)(q21q31)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = duplication_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_duplication_different_arms(self):
        """Duplication cannot have breakpoints on different arms."""
        from iscn_authenticator.rules.abnormality import duplication_breakpoint_rule

        bp1 = Breakpoint("p", 1, 2, None, False)
        bp2 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("dup", "1", [bp1, bp2], None, False, None, "dup(1)(p12q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = duplication_breakpoint_rule.validate(ast, abn)
        self.assertIn("same arm", errors[0].lower())

    def test_invalid_duplication_three_breakpoints(self):
        """Duplication cannot have three breakpoints."""
        from iscn_authenticator.rules.abnormality import duplication_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 2, 5, None, False)
        bp3 = Breakpoint("q", 3, 1, None, False)
        abn = Abnormality("dup", "1", [bp1, bp2, bp3], None, False, None, "dup(1)(q21q25q31)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = duplication_breakpoint_rule.validate(ast, abn)
        self.assertIn("one or two breakpoints", errors[0].lower())

    def test_skips_non_duplication(self):
        """Rule only applies to duplications."""
        from iscn_authenticator.rules.abnormality import duplication_breakpoint_rule

        bp1 = Breakpoint("q", 1, 3, None, False)
        bp2 = Breakpoint("q", 3, 3, None, False)
        abn = Abnormality("del", "5", [bp1, bp2], None, False, None, "del(5)(q13q33)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = duplication_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestFraBreakpointRule(unittest.TestCase):
    def test_valid_fra_one_breakpoint(self):
        """Fragile site with one breakpoint."""
        from iscn_authenticator.rules.abnormality import fra_breakpoint_rule

        bp1 = Breakpoint("q", 2, 7, 3, False)
        abn = Abnormality("fra", "X", [bp1], None, False, None, "fra(X)(q27.3)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = fra_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_fra_no_breakpoint(self):
        """Fragile site must have a breakpoint."""
        from iscn_authenticator.rules.abnormality import fra_breakpoint_rule

        abn = Abnormality("fra", "X", [], None, False, None, "fra(X)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = fra_breakpoint_rule.validate(ast, abn)
        self.assertIn("one breakpoint", errors[0].lower())

    def test_skips_non_fra(self):
        """Rule only applies to fragile sites."""
        from iscn_authenticator.rules.abnormality import fra_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("del", "5", [bp1], None, False, None, "del(5)(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = fra_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestInsBreakpointRule(unittest.TestCase):
    def test_valid_ins_three_breakpoints(self):
        """Insertion with three breakpoints."""
        from iscn_authenticator.rules.abnormality import ins_breakpoint_rule

        bp1 = Breakpoint("p", 1, 4, None, False)
        bp2 = Breakpoint("q", 2, 1, None, False)
        bp3 = Breakpoint("q", 3, 1, None, False)
        abn = Abnormality("ins", "5;2", [bp1, bp2, bp3], None, False, None, "ins(5;2)(p14;q21q31)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = ins_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_ins_two_breakpoints(self):
        """Insertion must have three breakpoints."""
        from iscn_authenticator.rules.abnormality import ins_breakpoint_rule

        bp1 = Breakpoint("p", 1, 4, None, False)
        bp2 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("ins", "5;2", [bp1, bp2], None, False, None, "ins(5;2)(p14;q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = ins_breakpoint_rule.validate(ast, abn)
        self.assertIn("three breakpoints", errors[0].lower())

    def test_skips_non_ins(self):
        """Rule only applies to insertions."""
        from iscn_authenticator.rules.abnormality import ins_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 3, 1, None, False)
        abn = Abnormality("dup", "1", [bp1, bp2], None, False, None, "dup(1)(q21q31)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = ins_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestDminBreakpointRule(unittest.TestCase):
    def test_valid_dmin_no_breakpoints(self):
        """Double minutes with no breakpoints."""
        from iscn_authenticator.rules.abnormality import dmin_breakpoint_rule

        abn = Abnormality("dmin", "", [], None, False, None, "dmin")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = dmin_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_dmin_with_breakpoint(self):
        """Double minutes should have no breakpoints."""
        from iscn_authenticator.rules.abnormality import dmin_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("dmin", "", [bp1], None, False, None, "dmin(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = dmin_breakpoint_rule.validate(ast, abn)
        self.assertIn("no breakpoints", errors[0].lower())

    def test_skips_non_dmin(self):
        """Rule only applies to double minutes."""
        from iscn_authenticator.rules.abnormality import dmin_breakpoint_rule

        abn = Abnormality("hsr", "", [], None, False, None, "hsr")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = dmin_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestHsrBreakpointRule(unittest.TestCase):
    def test_valid_hsr_no_breakpoints(self):
        """HSR with no breakpoints."""
        from iscn_authenticator.rules.abnormality import hsr_breakpoint_rule

        abn = Abnormality("hsr", "", [], None, False, None, "hsr")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = hsr_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_valid_hsr_one_breakpoint(self):
        """HSR with one breakpoint (location specified)."""
        from iscn_authenticator.rules.abnormality import hsr_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("hsr", "1", [bp1], None, False, None, "hsr(1)(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = hsr_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_hsr_two_breakpoints(self):
        """HSR should have at most one breakpoint."""
        from iscn_authenticator.rules.abnormality import hsr_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 3, 1, None, False)
        abn = Abnormality("hsr", "1", [bp1, bp2], None, False, None, "hsr(1)(q21q31)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = hsr_breakpoint_rule.validate(ast, abn)
        self.assertIn("zero or one breakpoint", errors[0].lower())

    def test_skips_non_hsr(self):
        """Rule only applies to HSR."""
        from iscn_authenticator.rules.abnormality import hsr_breakpoint_rule

        abn = Abnormality("dmin", "", [], None, False, None, "dmin")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = hsr_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestMarBreakpointRule(unittest.TestCase):
    def test_valid_mar_no_breakpoints(self):
        """Marker chromosome with no breakpoints."""
        from iscn_authenticator.rules.abnormality import mar_breakpoint_rule

        abn = Abnormality("mar", "", [], None, False, None, "+mar")
        ast = KaryotypeAST(47, "XX", [abn], None, None)
        errors = mar_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_mar_with_breakpoint(self):
        """Marker chromosome should have no breakpoints."""
        from iscn_authenticator.rules.abnormality import mar_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("mar", "", [bp1], None, False, None, "+mar(q21)")
        ast = KaryotypeAST(47, "XX", [abn], None, None)
        errors = mar_breakpoint_rule.validate(ast, abn)
        self.assertIn("no breakpoints", errors[0].lower())

    def test_skips_non_mar(self):
        """Rule only applies to marker chromosomes."""
        from iscn_authenticator.rules.abnormality import mar_breakpoint_rule

        abn = Abnormality("dmin", "", [], None, False, None, "dmin")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = mar_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestPseudicentricBreakpointRule(unittest.TestCase):
    def test_valid_pseudodicentric_two_chromosomes(self):
        """Pseudodicentric with two chromosomes and two breakpoints."""
        from iscn_authenticator.rules.abnormality import pseudodicentric_breakpoint_rule

        bp1 = Breakpoint("q", 1, 4, None, False)
        bp2 = Breakpoint("q", 1, 1, None, False)
        abn = Abnormality("psu dic", "13;14", [bp1, bp2], None, False, None, "psu dic(13;14)(q14;q11)")
        ast = KaryotypeAST(45, "XX", [abn], None, None)
        errors = pseudodicentric_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_pseudodicentric_mismatched(self):
        """Pseudodicentric with chromosome/breakpoint mismatch."""
        from iscn_authenticator.rules.abnormality import pseudodicentric_breakpoint_rule

        bp1 = Breakpoint("q", 1, 4, None, False)
        abn = Abnormality("psu dic", "13;14", [bp1], None, False, None, "psu dic(13;14)(q14)")
        ast = KaryotypeAST(45, "XX", [abn], None, None)
        errors = pseudodicentric_breakpoint_rule.validate(ast, abn)
        self.assertIn("2 chromosomes but 1 breakpoints", errors[0])

    def test_skips_non_pseudodicentric(self):
        """Rule only applies to pseudodicentric."""
        from iscn_authenticator.rules.abnormality import pseudodicentric_breakpoint_rule

        bp1 = Breakpoint("q", 1, 4, None, False)
        bp2 = Breakpoint("q", 1, 1, None, False)
        abn = Abnormality("dic", "13;14", [bp1, bp2], None, False, None, "dic(13;14)(q14;q11)")
        ast = KaryotypeAST(45, "XX", [abn], None, None)
        errors = pseudodicentric_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestAcentricBreakpointRule(unittest.TestCase):
    def test_valid_acentric_one_breakpoint(self):
        """Acentric with one breakpoint."""
        from iscn_authenticator.rules.abnormality import acentric_breakpoint_rule

        bp1 = Breakpoint("p", 1, 5, None, False)
        abn = Abnormality("ace", "5", [bp1], None, False, None, "ace(5)(p15)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = acentric_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_valid_acentric_two_breakpoints(self):
        """Acentric with two breakpoints."""
        from iscn_authenticator.rules.abnormality import acentric_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 3, 1, None, False)
        abn = Abnormality("ace", "1", [bp1, bp2], None, False, None, "ace(1)(q21q31)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = acentric_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_acentric_no_breakpoints(self):
        """Acentric must have at least one breakpoint."""
        from iscn_authenticator.rules.abnormality import acentric_breakpoint_rule

        abn = Abnormality("ace", "1", [], None, False, None, "ace(1)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = acentric_breakpoint_rule.validate(ast, abn)
        self.assertIn("1-2 breakpoints", errors[0])

    def test_invalid_acentric_three_breakpoints(self):
        """Acentric must have at most two breakpoints."""
        from iscn_authenticator.rules.abnormality import acentric_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 2, 5, None, False)
        bp3 = Breakpoint("q", 3, 1, None, False)
        abn = Abnormality("ace", "1", [bp1, bp2, bp3], None, False, None, "ace(1)(q21q25q31)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = acentric_breakpoint_rule.validate(ast, abn)
        self.assertIn("1-2 breakpoints", errors[0])

    def test_skips_non_acentric(self):
        """Rule only applies to acentric fragments."""
        from iscn_authenticator.rules.abnormality import acentric_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("del", "1", [bp1], None, False, None, "del(1)(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = acentric_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestTelomericAssociationBreakpointRule(unittest.TestCase):
    def test_valid_tas_two_chromosomes(self):
        """Telomeric association with two chromosomes and two breakpoints."""
        from iscn_authenticator.rules.abnormality import (
            telomeric_association_breakpoint_rule,
        )

        bp1 = Breakpoint("p", 1, 1, None, False)
        bp2 = Breakpoint("p", 1, 1, None, False)
        abn = Abnormality("tas", "13;14", [bp1, bp2], None, False, None, "tas(13;14)(p11;p11)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = telomeric_association_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_valid_tas_three_chromosomes(self):
        """Telomeric association with three chromosomes and three breakpoints."""
        from iscn_authenticator.rules.abnormality import (
            telomeric_association_breakpoint_rule,
        )

        bp1 = Breakpoint("p", 1, 1, None, False)
        bp2 = Breakpoint("p", 1, 1, None, False)
        bp3 = Breakpoint("p", 1, 1, None, False)
        abn = Abnormality(
            "tas",
            "13;14;15",
            [bp1, bp2, bp3],
            None,
            False,
            None,
            "tas(13;14;15)(p11;p11;p11)",
        )
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = telomeric_association_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_tas_mismatched(self):
        """Telomeric association with chromosome/breakpoint mismatch."""
        from iscn_authenticator.rules.abnormality import (
            telomeric_association_breakpoint_rule,
        )

        bp1 = Breakpoint("p", 1, 1, None, False)
        abn = Abnormality("tas", "13;14", [bp1], None, False, None, "tas(13;14)(p11)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = telomeric_association_breakpoint_rule.validate(ast, abn)
        self.assertIn("2 chromosomes but 1 breakpoints", errors[0])

    def test_skips_non_tas(self):
        """Rule only applies to telomeric associations."""
        from iscn_authenticator.rules.abnormality import (
            telomeric_association_breakpoint_rule,
        )

        bp1 = Breakpoint("q", 1, 4, None, False)
        bp2 = Breakpoint("q", 1, 1, None, False)
        abn = Abnormality("dic", "13;14", [bp1, bp2], None, False, None, "dic(13;14)(q14;q11)")
        ast = KaryotypeAST(45, "XX", [abn], None, None)
        errors = telomeric_association_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestFissionBreakpointRule(unittest.TestCase):
    def test_valid_fission_one_breakpoint(self):
        """Fission with one breakpoint."""
        from iscn_authenticator.rules.abnormality import fission_breakpoint_rule

        bp1 = Breakpoint("p", 1, 0, None, False)
        abn = Abnormality("fis", "1", [bp1], None, False, None, "fis(1)(p10)")
        ast = KaryotypeAST(47, "XX", [abn], None, None)
        errors = fission_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_fission_no_breakpoint(self):
        """Fission must have a breakpoint."""
        from iscn_authenticator.rules.abnormality import fission_breakpoint_rule

        abn = Abnormality("fis", "1", [], None, False, None, "fis(1)")
        ast = KaryotypeAST(47, "XX", [abn], None, None)
        errors = fission_breakpoint_rule.validate(ast, abn)
        self.assertIn("one breakpoint", errors[0].lower())

    def test_invalid_fission_two_breakpoints(self):
        """Fission should have only one breakpoint."""
        from iscn_authenticator.rules.abnormality import fission_breakpoint_rule

        bp1 = Breakpoint("p", 1, 0, None, False)
        bp2 = Breakpoint("q", 1, 0, None, False)
        abn = Abnormality("fis", "1", [bp1, bp2], None, False, None, "fis(1)(p10q10)")
        ast = KaryotypeAST(47, "XX", [abn], None, None)
        errors = fission_breakpoint_rule.validate(ast, abn)
        self.assertIn("one breakpoint", errors[0].lower())

    def test_skips_non_fission(self):
        """Rule only applies to fissions."""
        from iscn_authenticator.rules.abnormality import fission_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("del", "1", [bp1], None, False, None, "del(1)(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = fission_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestNeocentromereBreakpointRule(unittest.TestCase):
    def test_valid_neocentromere_one_breakpoint(self):
        """Neocentromere with one breakpoint."""
        from iscn_authenticator.rules.abnormality import neocentromere_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("neo", "1", [bp1], None, False, None, "neo(1)(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = neocentromere_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_neocentromere_no_breakpoint(self):
        """Neocentromere must have a breakpoint."""
        from iscn_authenticator.rules.abnormality import neocentromere_breakpoint_rule

        abn = Abnormality("neo", "1", [], None, False, None, "neo(1)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = neocentromere_breakpoint_rule.validate(ast, abn)
        self.assertIn("one breakpoint", errors[0].lower())

    def test_invalid_neocentromere_two_breakpoints(self):
        """Neocentromere should have only one breakpoint."""
        from iscn_authenticator.rules.abnormality import neocentromere_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        bp2 = Breakpoint("q", 2, 4, None, False)
        abn = Abnormality("neo", "1", [bp1, bp2], None, False, None, "neo(1)(q21q24)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = neocentromere_breakpoint_rule.validate(ast, abn)
        self.assertIn("one breakpoint", errors[0].lower())

    def test_skips_non_neocentromere(self):
        """Rule only applies to neocentromeres."""
        from iscn_authenticator.rules.abnormality import neocentromere_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("del", "1", [bp1], None, False, None, "del(1)(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = neocentromere_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


class TestIncompleteBreakpointRule(unittest.TestCase):
    def test_valid_incomplete_no_breakpoints(self):
        """Incomplete marker with no breakpoints."""
        from iscn_authenticator.rules.abnormality import incomplete_breakpoint_rule

        abn = Abnormality("inc", "", [], None, False, None, "inc")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = incomplete_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])

    def test_invalid_incomplete_with_breakpoint(self):
        """Incomplete marker should have no breakpoints."""
        from iscn_authenticator.rules.abnormality import incomplete_breakpoint_rule

        bp1 = Breakpoint("q", 2, 1, None, False)
        abn = Abnormality("inc", "", [bp1], None, False, None, "inc(q21)")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = incomplete_breakpoint_rule.validate(ast, abn)
        self.assertIn("no breakpoints", errors[0].lower())

    def test_skips_non_incomplete(self):
        """Rule only applies to incomplete markers."""
        from iscn_authenticator.rules.abnormality import incomplete_breakpoint_rule

        abn = Abnormality("dmin", "", [], None, False, None, "dmin")
        ast = KaryotypeAST(46, "XX", [abn], None, None)
        errors = incomplete_breakpoint_rule.validate(ast, abn)
        self.assertEqual(errors, [])


if __name__ == "__main__":
    unittest.main()
