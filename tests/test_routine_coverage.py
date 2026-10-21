import json
import unittest
from pathlib import Path

from iscn_authenticator import validate_karyotype

CASES = json.loads(Path(__file__).with_name("routine_cases.json").read_text())


class RoutineCoverageTests(unittest.TestCase):
    def test_accepted(self):
        for notation in CASES["valid"]:
            with self.subTest(notation=notation):
                result = validate_karyotype(notation)
                self.assertTrue(result.valid, result.errors)

    def test_rejected(self):
        for notation in CASES["invalid"]:
            with self.subTest(notation=notation):
                self.assertFalse(validate_karyotype(notation).valid)

    def test_coverage_metadata(self):
        import iscn_authenticator

        coverage = getattr(iscn_authenticator, "get_validation_coverage", None)
        self.assertTrue(callable(coverage))
        self.assertEqual(coverage(validate_karyotype("46,XX[20]").parsed)["level"], "supported")
        for text in [
            "46,XX,der(2)t(2;8)(p21;q24)",
            "46,XX,foo",
            "69,XXY",
            "45~48,XX,+7[cp10]",
            "46,XX.arr",
            "nuc ish (ABC)x2[20]",
            "46,XX,del(5)(p15.311)",
            "46,XX,i(X)(q10)",
            "46,XX,dmin",
        ]:
            with self.subTest(notation=text):
                self.assertEqual(coverage(validate_karyotype(text).parsed)["level"], "partial")

    def test_raw_count_and_clone_provenance(self):
        single = validate_karyotype(" 47,XX,+21[20] ").parsed
        self.assertEqual(single.raw, " 47,XX,+21[20] ")
        self.assertEqual(single.count, 20)
        self.assertIsNone(single.cell_lines)
        raw = "46,XX,t(2;8)(p21;q24)[7]/47,sl,+12[4]/48,sdl1,+9[3]"
        line = validate_karyotype(raw).parsed.cell_lines[2]
        self.assertEqual(line.raw, "48,sdl1,+9[3]")
        self.assertEqual(line.reference.token, "sdl1")
        self.assertEqual(line.reference.cell_line, 2)
        self.assertEqual([a.raw for a in line.changes], ["+9"])
        self.assertEqual([a.raw for a in line.abnormalities], ["t(2;8)(p21;q24)", "+12", "+9"])

    def test_oversized_numbers_do_not_crash(self):
        for text in ["46,XX,del(5)(p15)[2]/47,sdl" + "9" * 4500 + ",+8[1]", "46,XX,del(5)(p15)×" + "9" * 4500]:
            self.assertFalse(validate_karyotype(text).valid)
        self.assertTrue(validate_karyotype("0" * 4500 + "45~46,XX").valid)

    def test_region_is_not_centromere(self):
        result = validate_karyotype("46,XX,del(5)(p1?)")
        self.assertTrue(result.valid)
        self.assertIsNone(result.parsed.abnormalities[0].breakpoints[0].band)
        self.assertTrue(result.parsed.abnormalities[0].breakpoints[0].uncertain)
