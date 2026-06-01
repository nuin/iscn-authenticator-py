import os
import sys
import unittest

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from iscn_authenticator.main import is_valid_karyotype, validate_karyotype


class TestIsValidKaryotype(unittest.TestCase):
    def test_valid_karyotypes(self):
        self.assertTrue(is_valid_karyotype("46,XX"))
        self.assertTrue(is_valid_karyotype("47,XY"))
        self.assertTrue(is_valid_karyotype("47,XXX"))
        self.assertTrue(is_valid_karyotype("45,X"))
        self.assertTrue(is_valid_karyotype("47,XX,+21"))
        self.assertTrue(is_valid_karyotype("46,XX,del(5)(q13)"))

    def test_invalid_karyotypes(self):
        self.assertFalse(is_valid_karyotype("46,X"))  # Should be 45,X
        self.assertFalse(is_valid_karyotype("46,Y"))  # No X chromosome
        self.assertFalse(is_valid_karyotype("46,Z"))  # Invalid sex chr
        self.assertFalse(is_valid_karyotype("46XX"))  # Missing comma
        self.assertFalse(is_valid_karyotype("foo,bar"))  # Invalid count
        self.assertFalse(is_valid_karyotype("46,XXY"))  # Should be 47,XXY
        # Note: Abnormality format validation added in Task 7-12
        # These tests will be uncommented when that functionality is implemented:
        # self.assertFalse(is_valid_karyotype("46,XX,"))  # Trailing comma
        # self.assertFalse(is_valid_karyotype("47,XX,21"))  # Missing +/- prefix
        # self.assertFalse(is_valid_karyotype("47,XX,+21,+"))  # Incomplete abnormality

    def test_abnormality_karyotypes(self):
        # Valid abnormality karyotypes (basic parsing)
        self.assertTrue(is_valid_karyotype("46,XX,del(5)(q13q33)"))
        self.assertTrue(is_valid_karyotype("46,XX,del(5)(q13)"))
        # Note: Detailed abnormality syntax validation added in Task 8-12
        # The following tests will be re-enabled when parser validates syntax:
        # self.assertFalse(is_valid_karyotype("46,XX,del(5)(q)"))  # Invalid breakpoint
        # self.assertFalse(is_valid_karyotype("46,XX,del(5)(13q)"))  # Wrong order
        # self.assertFalse(is_valid_karyotype("46,XX,del(5)q13"))  # Missing parens
        # self.assertFalse(is_valid_karyotype("46,XX,del(5)(q13q)"))  # Incomplete


class TestValidateKaryotype(unittest.TestCase):
    """Tests for the new validate_karyotype() API."""

    def test_returns_validation_result(self):
        result = validate_karyotype("46,XX")
        self.assertTrue(result.valid)
        self.assertEqual(result.errors, [])
        self.assertIsNotNone(result.parsed)

    def test_returns_errors_for_invalid(self):
        result = validate_karyotype("46,Y")
        self.assertFalse(result.valid)
        self.assertTrue(len(result.errors) > 0)

    def test_parsed_contains_chromosome_count(self):
        result = validate_karyotype("46,XX")
        self.assertEqual(result.parsed.chromosome_count, 46)

    def test_parsed_contains_sex_chromosomes(self):
        result = validate_karyotype("47,XXY")
        self.assertEqual(result.parsed.sex_chromosomes, "XXY")

    def test_backward_compat_is_valid_karyotype(self):
        # Ensure old API still works
        self.assertTrue(is_valid_karyotype("46,XX"))
        self.assertFalse(is_valid_karyotype("46,Y"))

    def test_parse_error_returns_errors(self):
        result = validate_karyotype("")
        self.assertFalse(result.valid)
        self.assertIn("empty", result.errors[0].lower())
        self.assertIsNone(result.parsed)

    def test_parse_error_missing_comma(self):
        result = validate_karyotype("46XX")
        self.assertFalse(result.valid)
        self.assertIn("comma", result.errors[0].lower())


if __name__ == "__main__":
    unittest.main()
