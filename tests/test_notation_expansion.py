"""Synthetic regression controls for bounded complex notation checks."""

import unittest

from iscn_authenticator.main import validate_karyotype


class TestNotationExpansion(unittest.TestCase):
    def test_nested_derivative_components(self):
        cases = [
            ("46,XX,der(5)t(5;99)(q13;p13)", False, "Invalid chromosome '99'"),
            ("46,XX,der(5)t(5;12)(q13;p99)", False, "12p99"),
            ("46,XX,der(5)t(5;12)(q13)", False, "breakpoints"),
            ("46,XX,der(5)t(5;12)(q13;p13", False, "parentheses"),
            ("46,XX,der(5)t(5;12)(q13;p13)", True, None),
            ("47,XX,+der(5)t(5;12)(q13;p13)", True, None),
            ("46,XX,der(5)t(5;99)", False, "Invalid chromosome '99'"),
            ("46,XX,der(5)t(5;12)", True, None),
            ("46,XX,der(11)t(5;11)(q13;q12)t(5;8)(q23;q24.1)", True, None),
            ("45,XX,der(21;21)(q10;q10)", True, None),
            ("45,XX,der(21;21)(q10;q99)", False, "21q99"),
            ("45,XX,der(21;21)(q10)", False, "breakpoints"),
        ]
        for notation, valid, error in cases:
            with self.subTest(notation=notation):
                result = validate_karyotype(notation)
                self.assertEqual(result.valid, valid, result.errors)
                if error:
                    self.assertIn(error, "; ".join(result.errors))
                self.assertEqual(result.parsed.raw, notation)

    def test_composite_counts_are_group_local(self):
        result = validate_karyotype("45~47,XX,del(5)(q13q33),+8[cp12]/46,XX[8]")
        self.assertTrue(result.valid, result.errors)
        self.assertEqual([line.count for line in result.parsed.cell_lines], [12, 8])
        self.assertIsNone(result.parsed.modifiers)
        result = validate_karyotype("45~47,XX,+8[cp12]/47,XX,+8,+9[8]")
        self.assertFalse(result.valid)
        self.assertIn("Cell line 2: Chromosome count 47 does not match the expected total 48", "; ".join(result.errors))

    def test_single_composite_count(self):
        for suffix in ("[cp12]", ".cp[12]"):
            result = validate_karyotype(f"45~47,XX,+8{suffix}")
            self.assertTrue(result.valid, result.errors)
            self.assertEqual(result.parsed.count, 12)
            self.assertEqual(result.parsed.modifiers.composite, 12)

    def test_composite_reference_rejected(self):
        for notation in ("47,XX,+8[cp12]/48,idem,+9[8]", "47,XX,+8[12]/48,idem,+9[cp8]"):
            result = validate_karyotype(notation)
            self.assertFalse(result.valid)
            self.assertIn("Composite groups cannot use clone references", "; ".join(result.errors))

    def test_clone_arithmetic(self):
        for notation, sex, findings in (
            ("47,XY,+8[10]/94,idemx2[5]", "XXYY", [("+", 2)]),
            ("46,XX,inv(3)(p21q26)[10]/45,idem,-inv(3)[5]", "XX", [("-", None)]),
            ("46,XX,inv(3)(p21q26)[10]/46,idem,+3,-inv(3)[5]", "XX", []),
            ("47,XX,+8[10]/46,idem,-8[5]", "XX", []),
        ):
            with self.subTest(notation=notation):
                result = validate_karyotype(notation)
                self.assertTrue(result.valid, result.errors)
                line = result.parsed.cell_lines[1]
                self.assertEqual(line.sex_chromosomes, sex)
                self.assertEqual([(abn.type, abn.copy_count) for abn in line.abnormalities], findings)
                self.assertIsNone(result.parsed.cell_lines[0].abnormalities[0].copy_count)

    def test_unsafe_clone_arithmetic(self):
        for notation, error in (
            (
                "46,XX,inv(3)(p21q26)[5]/93,idemx2,+8x9007199254740991,+9x2,-8x9007199254740991[5]",
                "Clone count change exceeds the safe integer limit",
            ),
            ("47,XX,+8[5]/93,idemx2[5]", "expected total 94"),
            ("46,XX,inv(3)(p21q26)[5]/46,idem,-inv(3)[5]", "expected total 45"),
            ("46,XX,inv(3)(p21q26),inv(3)(p13q21)[5]/45,idem,-inv(3)[5]", "ambiguous"),
            ("46,XX,inv(3)(p21q26)[5]/45,idem,-inv(5)[5]", "absent"),
            ("46,XX,?inv(3)(p21q26)[5]/45,idem,-inv(3)[5]", "uncertain"),
            ("47,XX,+8[5]/94,idemx0[5]", "multiplier"),
            ("47,XX,+8[5]/94,idemx9007199254740992[5]", "multiplier"),
            ("47,XX,+8[5]/94,idemx1000000000[5]", "expansion limit"),
        ):
            with self.subTest(notation=notation):
                result = validate_karyotype(notation)
                self.assertFalse(result.valid)
                self.assertIn(error, "; ".join(result.errors))
