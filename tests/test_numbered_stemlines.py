"""Reference-selection regressions based on ISCN 2024 section 6.3.4."""

import unittest

from iscn_authenticator import validate_karyotype


class TestNumberedStemlines(unittest.TestCase):
    def test_separate_families_and_multipliers(self):
        notation = "45,XY,-7[5]/46,sl1,+8[6]/46,XY,t(9;22)(q34;q11.2)[3]/92,sl2×2[5]/93,sl2x2,+8[2]"
        result = validate_karyotype(notation)
        self.assertTrue(result.valid, result.errors)
        lines = result.parsed.cell_lines
        sources = [line.reference.cell_line if line.reference else None for line in lines]
        self.assertEqual(sources, [None, 1, None, 3, 3])
        self.assertEqual([line.sex_chromosomes for line in lines], ["XY", "XY", "XY", "XXYY", "XXYY"])
        self.assertEqual([[(a.type, a.chromosome, a.copy_count) for a in line.abnormalities] for line in lines], [
            [("-", "7", None)], [("-", "7", None), ("+", "8", None)], [("t", "9;22", None)],
            [("t", "9;22", 2)], [("t", "9;22", 2), ("+", "8", None)],
        ])
        self.assertEqual(lines[4].raw, "93,sl2x2,+8[2]")
        self.assertEqual(lines[4].reference.token, "sl2x2")
        self.assertEqual([a.raw for a in lines[4].changes], ["+8"])

    def test_references_select_the_intended_family(self):
        for notation, index, sex, findings in (
            ("47,XY,+8[5]/47,XX,+9[5]/48,sl2,+10[5]/48,sl1,+11[5]", 1, "XY", ["+8", "+11"]),
            ("47,XY,+8[5]/47,XX,+9[5]/48,idem,+10[5]", 1, "XY", ["+8", "+10"]),
            ("47,XX,+8[5]//47,XY,+9[5]/47,XY,+10[5]/48,sl2,+11[5]", 2, "XY", ["+10", "+11"]),
            ("47,XX,+8[cp5]/47,XY,+9[5]/48,sl2,+10[5]", 2, "XY", ["+9", "+10"]),
        ):
            with self.subTest(notation=notation):
                result = validate_karyotype(notation)
                self.assertTrue(result.valid, result.errors)
                line = result.parsed.cell_lines[-1]
                self.assertEqual(line.reference.cell_line, index)
                self.assertEqual(line.sex_chromosomes, sex)
                self.assertEqual([a.raw for a in line.abnormalities], findings)
                self.assertEqual(line.is_donor, "//" in notation)

    def test_invalid_or_ambiguous_selection(self):
        for notation, error in (
            ("47,XX,+8[5]/48,sl,+9[5]/47,XX,+10[5]/49,sdl1,+11[5]", "Sideline references with multiple stemlines"),
            ("47,XX,+8[5]/47,XX,+9[5]/48,sl,+10[5]", "ambiguous after an independent clone"),
            ("47,XX,+8[5]/48,sl2,+9[5]/47,XX,+10[5]", "no preceding stemline"),
            ("47,XX,+8[5]/48,sl0,+9[5]", "Unsupported clone reference"),
            ("47,XX,+8[5]/48,sl9007199254740992,+9[5]", "no preceding stemline"),
            ("47,XX,+8[5]/48,sl1x0,+9[5]", "multiplier"),
            ("47,XX,+8[5]//48,sl1,+9[5]", "no preceding reference clone"),
            ("47,XX,+8[cp5]/47,XX,+9[5]/48,sl1,+10[5]", "Composite groups cannot use clone references"),
            ("47,XX,+8[5]/48,sl1,+9[5]/49,idem,+10[5]", "Do not mix"),
            ("47,XX,+8[5]/93,sl1x2[5]", "expected total 94"),
            ("47,XX,+8[5]/48,sl" + "9" * 5000 + ",+9[5]", "no preceding stemline"),
        ):
            with self.subTest(notation=notation[:200]):
                result = validate_karyotype(notation)
                self.assertFalse(result.valid)
                self.assertIn(error, "; ".join(result.errors))
