import unittest

from iscn_authenticator.explain import explain, generate_template_explanation
from iscn_authenticator.parser import KaryotypeParser


class TestExplain(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parser = KaryotypeParser()

    def test_template_normal_female(self):
        ast = self.parser.parse("46,XX")
        result = generate_template_explanation(ast)
        self.assertEqual(result.summary, "46,XX karyotype with 0 abnormalities.")
        self.assertIn("total chromosome count of 46", result.detail)
        self.assertEqual(result.confidence, "template")

    def test_template_trisomy_21(self):
        ast = self.parser.parse("47,XY,+21")
        abn = ast.abnormalities[0]
        result = generate_template_explanation(abn)
        self.assertEqual(result.summary, "Gain of chromosome 21.")
        self.assertEqual(result.confidence, "template")

    def test_template_deletion_with_breakpoints(self):
        ast = self.parser.parse("46,XX,del(5)(q13q33)")
        abn = ast.abnormalities[0]
        result = generate_template_explanation(abn)
        self.assertEqual(result.summary, "Deletion on chromosome 5 at q13, q33.")
        self.assertEqual(result.confidence, "template")

    def test_explain_returns_template(self):
        # v0.2.0 always returns a template explanation; curated lookup
        # is deferred to a future release.
        ast = self.parser.parse("46,XX,del(7)(q22q36)")
        abn = ast.abnormalities[0]
        result = explain(abn)
        self.assertEqual(result.confidence, "template")
        self.assertEqual(result.summary, "Deletion on chromosome 7 at q22, q36.")


if __name__ == "__main__":
    unittest.main()
