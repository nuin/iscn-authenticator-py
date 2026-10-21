"""Synthetic controls for interphase prefix routing, not complete FISH grammar."""

import unittest

from iscn_authenticator import get_validation_coverage, validate_karyotype


class TestInterphasePrefix(unittest.TestCase):
    def test_spacing_preserves_payload_without_creating_cell_lines(self):
        payload = "(D21S259)x3[25/100]"
        for separator in ("", " ", "\t"):
            with self.subTest(separator=separator):
                raw = f"nuc ish{separator}{payload}"
                result = validate_karyotype(raw)
                self.assertTrue(result.valid, result.errors)
                self.assertEqual(result.parsed.raw, raw)
                self.assertEqual(result.parsed.modifiers.ish, payload)
                self.assertTrue(result.parsed.modifiers.interphase)
                self.assertIsNone(result.parsed.cell_lines)
                self.assertIsNone(result.parsed.chromosome_count)
                self.assertEqual(get_validation_coverage(result.parsed)["level"], "partial")

    def test_similar_but_invalid_prefix_is_not_accepted(self):
        self.assertFalse(validate_karyotype("nuc ishx(D21S259)x3[25/100]").valid)
