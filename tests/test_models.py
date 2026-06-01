# tests/test_models.py
import unittest

from iscn_authenticator.models import (
    Abnormality,
    Breakpoint,
    KaryotypeAST,
    ValidationResult,
)


class TestValidationResult(unittest.TestCase):
    def test_valid_result(self):
        result = ValidationResult(valid=True, errors=[], parsed=None)
        self.assertTrue(result.valid)
        self.assertEqual(result.errors, [])

    def test_invalid_result_with_errors(self):
        result = ValidationResult(valid=False, errors=["Invalid chromosome count"], parsed=None)
        self.assertFalse(result.valid)
        self.assertEqual(result.errors, ["Invalid chromosome count"])


class TestBreakpoint(unittest.TestCase):
    def test_breakpoint_basic(self):
        bp = Breakpoint(arm="q", region=13, band=3, subband=None, uncertain=False)
        self.assertEqual(bp.arm, "q")
        self.assertEqual(bp.region, 13)
        self.assertEqual(bp.band, 3)

    def test_breakpoint_with_subband(self):
        bp = Breakpoint(arm="p", region=11, band=2, subband="1", uncertain=False)
        self.assertEqual(bp.subband, "1")


class TestAbnormality(unittest.TestCase):
    def test_numerical_abnormality(self):
        abn = Abnormality(
            type="+",
            chromosome="21",
            breakpoints=[],
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw="+21",
        )
        self.assertEqual(abn.type, "+")
        self.assertEqual(abn.chromosome, "21")


class TestKaryotypeAST(unittest.TestCase):
    def test_simple_karyotype(self):
        ast = KaryotypeAST(
            chromosome_count=46,
            sex_chromosomes="XX",
            abnormalities=[],
            cell_lines=None,
            modifiers=None,
        )
        self.assertEqual(ast.chromosome_count, 46)
        self.assertEqual(ast.sex_chromosomes, "XX")


if __name__ == "__main__":
    unittest.main()
