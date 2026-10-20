import copy
import unittest

from iscn_authenticator.engine import RuleEngine
from iscn_authenticator.main import validate_karyotype
from iscn_authenticator.parser import KaryotypeParser
from iscn_authenticator.rules.base import Rule


class TestCellLines(unittest.TestCase):
    def test_later_line_errors_have_positions(self):
        cases = [
            ("47,XX,+21[10]/47,XX,+99[5]", "Cell line 2: Invalid chromosome '99'"),
            ("47,XX,+21[5]/46,Y[10]", "Cell line 2: Sex chromosomes"),
            ("47,XX,+21[5]/45,X[5]/0,XX[10]", "Cell line 3: Chromosome count 0"),
        ]
        for notation, error in cases:
            with self.subTest(notation=notation):
                result = validate_karyotype(notation)
                self.assertFalse(result.valid)
                self.assertIn(error, "\n".join(result.errors))

    def test_engine_runs_once_per_line_without_mutating_ast(self):
        ast = KaryotypeParser().parse("47,XX,+21[5]/46,XX[10]")
        before = copy.deepcopy(ast)
        engine = RuleEngine()
        engine.add_rule(
            Rule(
                id="context",
                category="test",
                description="Current cell line",
                validate=lambda line, target: [f"{line.chromosome_count},{line.sex_chromosomes}"],
            )
        )
        result = engine.validate(ast)
        self.assertEqual(result.errors, ["Cell line 1: 47,XX", "Cell line 2: 46,XX"])
        self.assertEqual(ast, before)
        self.assertIs(result.parsed, ast)

    def test_prefix_range_and_donor_information_survive_parsing(self):
        parser = KaryotypeParser()
        self.assertTrue(parser.parse("mos 47,XX,+21[5]/46,XX[10]").modifiers.mosaic)
        self.assertEqual(parser.parse("47,XX,+21[5]/45~48,XX[10]").cell_lines[1].chromosome_count, "45~48")
        ast = parser.parse("46,XX[5]//47,XY,+21[15]/46,XY[2]")
        self.assertEqual([line.is_donor for line in ast.cell_lines], [False, True, True])

    def test_missing_and_malformed_counts_have_positions(self):
        for suffix in ["", "[0]", "[-1]", "[1.5]", "[١٠]", "[9007199254740992]", "[" + "9" * 5000 + "]"]:
            with self.subTest(suffix=suffix):
                result = validate_karyotype(f"47,XX,+21[5]/46,XX{suffix}")
                self.assertFalse(result.valid)
                self.assertIn("Cell line 2: Cell count", "\n".join(result.errors))
