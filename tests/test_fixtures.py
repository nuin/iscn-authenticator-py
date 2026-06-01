"""Fixture-driven cross-implementation consistency tests.

Reads fixtures/validity.json and asserts that iscn_authenticator agrees with
each case. A parallel TypeScript runner (packages/core/tests/fixtures.test.ts)
reads the same file. CI runs both.
"""

import json
import pathlib
import unittest

from iscn_authenticator.main import is_valid_karyotype

FIXTURES_PATH = pathlib.Path(__file__).resolve().parent.parent / "fixtures" / "validity.json"


class TestValidityFixtures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with FIXTURES_PATH.open() as f:
            cls.fixtures = json.load(f)

    def test_valid_cases_all_pass(self):
        for case in self.fixtures["valid"]:
            with self.subTest(input=case["input"], note=case.get("note", "")):
                self.assertTrue(
                    is_valid_karyotype(case["input"]),
                    f"expected valid: {case['input']!r} ({case.get('note', '')})",
                )

    def test_invalid_cases_all_fail(self):
        for case in self.fixtures["invalid"]:
            with self.subTest(input=case["input"], reason=case.get("reason", "")):
                self.assertFalse(
                    is_valid_karyotype(case["input"]),
                    f"expected invalid: {case['input']!r} ({case.get('reason', '')})",
                )


if __name__ == "__main__":
    unittest.main()
