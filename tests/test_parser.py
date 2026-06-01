# tests/test_parser.py
import unittest

from iscn_authenticator.parser import KaryotypeParser, ParseError


class TestKaryotypeParserBasic(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_normal_female(self):
        result = self.parser.parse("46,XX")
        self.assertEqual(result.chromosome_count, 46)
        self.assertEqual(result.sex_chromosomes, "XX")
        self.assertEqual(result.abnormalities, [])

    def test_parse_normal_male(self):
        result = self.parser.parse("46,XY")
        self.assertEqual(result.chromosome_count, 46)
        self.assertEqual(result.sex_chromosomes, "XY")

    def test_parse_turner_syndrome(self):
        result = self.parser.parse("45,X")
        self.assertEqual(result.chromosome_count, 45)
        self.assertEqual(result.sex_chromosomes, "X")

    def test_parse_klinefelter(self):
        result = self.parser.parse("47,XXY")
        self.assertEqual(result.chromosome_count, 47)
        self.assertEqual(result.sex_chromosomes, "XXY")

    def test_parse_triple_x(self):
        result = self.parser.parse("47,XXX")
        self.assertEqual(result.chromosome_count, 47)
        self.assertEqual(result.sex_chromosomes, "XXX")

    def test_parse_xyy(self):
        result = self.parser.parse("47,XYY")
        self.assertEqual(result.chromosome_count, 47)
        self.assertEqual(result.sex_chromosomes, "XYY")

    def test_parse_undisclosed_sex(self):
        result = self.parser.parse("46,U")
        self.assertEqual(result.chromosome_count, 46)
        self.assertEqual(result.sex_chromosomes, "U")

    def test_parse_empty_string_raises(self):
        with self.assertRaises(ParseError) as ctx:
            self.parser.parse("")
        self.assertIn("empty", str(ctx.exception).lower())

    def test_parse_missing_comma_raises(self):
        with self.assertRaises(ParseError) as ctx:
            self.parser.parse("46XX")
        self.assertIn("comma", str(ctx.exception).lower())

    def test_parse_invalid_count_raises(self):
        with self.assertRaises(ParseError) as ctx:
            self.parser.parse("foo,XX")
        self.assertIn("chromosome count", str(ctx.exception).lower())

    def test_parse_whitespace_handling(self):
        result = self.parser.parse("  46 , XX  ")
        self.assertEqual(result.chromosome_count, 46)
        self.assertEqual(result.sex_chromosomes, "XX")

    def test_parse_range_notation(self):
        result = self.parser.parse("45~48,XX")
        self.assertEqual(result.chromosome_count, "45~48")
        self.assertEqual(result.sex_chromosomes, "XX")


class TestParserNumericalAbnormalities(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_trisomy_21(self):
        result = self.parser.parse("47,XX,+21")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "+")
        self.assertEqual(abn.chromosome, "21")

    def test_parse_monosomy_7(self):
        result = self.parser.parse("45,XY,-7")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "-")
        self.assertEqual(abn.chromosome, "7")

    def test_parse_multiple_numerical(self):
        result = self.parser.parse("48,XY,+18,+21")
        self.assertEqual(len(result.abnormalities), 2)
        self.assertEqual(result.abnormalities[0].chromosome, "18")
        self.assertEqual(result.abnormalities[1].chromosome, "21")

    def test_parse_sex_chromosome_gain(self):
        result = self.parser.parse("48,XXXY,+X")
        self.assertEqual(len(result.abnormalities), 1)
        self.assertEqual(result.abnormalities[0].chromosome, "X")

    def test_parse_sex_chromosome_loss(self):
        result = self.parser.parse("45,XY,-Y")
        self.assertEqual(len(result.abnormalities), 1)
        self.assertEqual(result.abnormalities[0].type, "-")
        self.assertEqual(result.abnormalities[0].chromosome, "Y")


class TestParserDeletions(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_terminal_deletion(self):
        result = self.parser.parse("46,XX,del(5)(q13)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "del")
        self.assertEqual(abn.chromosome, "5")
        self.assertEqual(len(abn.breakpoints), 1)
        self.assertEqual(abn.breakpoints[0].arm, "q")
        self.assertEqual(abn.breakpoints[0].region, 1)
        self.assertEqual(abn.breakpoints[0].band, 3)

    def test_parse_interstitial_deletion(self):
        result = self.parser.parse("46,XX,del(5)(q13q33)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "del")
        self.assertEqual(len(abn.breakpoints), 2)
        self.assertEqual(abn.breakpoints[0].arm, "q")
        self.assertEqual(abn.breakpoints[1].arm, "q")
        self.assertEqual(abn.breakpoints[1].region, 3)
        self.assertEqual(abn.breakpoints[1].band, 3)

    def test_parse_deletion_with_subband(self):
        result = self.parser.parse("46,XY,del(7)(p11.2)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.breakpoints[0].subband, "2")
        self.assertEqual(abn.breakpoints[0].region, 1)
        self.assertEqual(abn.breakpoints[0].band, 1)

    def test_parse_deletion_x_chromosome(self):
        result = self.parser.parse("46,X,del(X)(p22)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.chromosome, "X")
        self.assertEqual(abn.breakpoints[0].arm, "p")
        self.assertEqual(abn.breakpoints[0].region, 2)
        self.assertEqual(abn.breakpoints[0].band, 2)

    def test_parse_deletion_p_arm(self):
        result = self.parser.parse("46,XY,del(1)(p36)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.breakpoints[0].arm, "p")
        self.assertEqual(abn.breakpoints[0].region, 3)
        self.assertEqual(abn.breakpoints[0].band, 6)


class TestParserDuplications(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_duplication(self):
        result = self.parser.parse("46,XX,dup(1)(p31p22)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "dup")
        self.assertEqual(abn.chromosome, "1")
        self.assertEqual(len(abn.breakpoints), 2)

    def test_parse_tandem_duplication(self):
        result = self.parser.parse("46,XY,dup(7)(q11.2q22)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.breakpoints[0].arm, "q")
        self.assertEqual(abn.breakpoints[1].arm, "q")
        self.assertEqual(abn.breakpoints[0].subband, "2")

    def test_parse_single_breakpoint_duplication(self):
        result = self.parser.parse("46,XX,dup(3)(q21)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "dup")
        self.assertEqual(len(abn.breakpoints), 1)
        self.assertEqual(abn.breakpoints[0].arm, "q")


class TestParserInversions(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_pericentric_inversion(self):
        result = self.parser.parse("46,XX,inv(9)(p12q13)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "inv")
        self.assertEqual(abn.chromosome, "9")
        self.assertEqual(len(abn.breakpoints), 2)
        self.assertEqual(abn.breakpoints[0].arm, "p")
        self.assertEqual(abn.breakpoints[1].arm, "q")

    def test_parse_paracentric_inversion(self):
        result = self.parser.parse("46,XY,inv(3)(q21q26)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.breakpoints[0].arm, "q")
        self.assertEqual(abn.breakpoints[1].arm, "q")
        self.assertEqual(abn.breakpoints[0].region, 2)
        self.assertEqual(abn.breakpoints[0].band, 1)
        self.assertEqual(abn.breakpoints[1].region, 2)
        self.assertEqual(abn.breakpoints[1].band, 6)

    def test_parse_inversion_with_subband(self):
        result = self.parser.parse("46,XY,inv(2)(p11.2q13)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.breakpoints[0].subband, "2")


class TestParserTranslocations(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_reciprocal_translocation(self):
        result = self.parser.parse("46,XX,t(9;22)(q34;q11.2)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "t")
        self.assertEqual(abn.chromosome, "9;22")
        self.assertEqual(len(abn.breakpoints), 2)
        self.assertEqual(abn.breakpoints[0].arm, "q")
        self.assertEqual(abn.breakpoints[1].arm, "q")
        self.assertEqual(abn.breakpoints[1].subband, "2")

    def test_parse_three_way_translocation(self):
        result = self.parser.parse("46,XY,t(1;3;5)(p32;q21;q31)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.chromosome, "1;3;5")
        self.assertEqual(len(abn.breakpoints), 3)
        self.assertEqual(abn.breakpoints[0].arm, "p")
        self.assertEqual(abn.breakpoints[1].arm, "q")
        self.assertEqual(abn.breakpoints[2].arm, "q")

    def test_parse_translocation_sex_chromosome(self):
        result = self.parser.parse("46,X,t(X;18)(p11.2;q21)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.chromosome, "X;18")
        self.assertEqual(len(abn.breakpoints), 2)


class TestParserIsochromosomes(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_isochromosome_short_form(self):
        """Test i(17q) - short form with arm in parentheses."""
        result = self.parser.parse("46,XX,i(17q)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "i")
        self.assertEqual(abn.chromosome, "17")
        self.assertEqual(len(abn.breakpoints), 1)
        self.assertEqual(abn.breakpoints[0].arm, "q")

    def test_parse_isochromosome_long_form(self):
        """Test i(17)(q10) - long form with breakpoint."""
        result = self.parser.parse("46,XX,i(17)(q10)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "i")
        self.assertEqual(abn.chromosome, "17")
        self.assertEqual(len(abn.breakpoints), 1)
        self.assertEqual(abn.breakpoints[0].arm, "q")
        self.assertEqual(abn.breakpoints[0].region, 1)
        self.assertEqual(abn.breakpoints[0].band, 0)

    def test_parse_isochromosome_p_arm(self):
        """Test i(9p) - isochromosome of p arm."""
        result = self.parser.parse("46,XY,i(9p)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "i")
        self.assertEqual(abn.chromosome, "9")
        self.assertEqual(abn.breakpoints[0].arm, "p")

    def test_parse_isochromosome_x_chromosome(self):
        """Test i(Xq) - isochromosome of X chromosome q arm."""
        result = self.parser.parse("46,X,i(Xq)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "i")
        self.assertEqual(abn.chromosome, "X")
        self.assertEqual(abn.breakpoints[0].arm, "q")


class TestParserRingChromosomes(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_ring_simple(self):
        """Test r(1) - simple ring chromosome."""
        result = self.parser.parse("46,XX,r(1)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "r")
        self.assertEqual(abn.chromosome, "1")
        self.assertEqual(len(abn.breakpoints), 0)

    def test_parse_ring_with_breakpoints(self):
        """Test r(1)(p36q42) - ring with breakpoints."""
        result = self.parser.parse("46,XY,r(1)(p36q42)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "r")
        self.assertEqual(abn.chromosome, "1")
        self.assertEqual(len(abn.breakpoints), 2)
        self.assertEqual(abn.breakpoints[0].arm, "p")
        self.assertEqual(abn.breakpoints[0].region, 3)
        self.assertEqual(abn.breakpoints[0].band, 6)
        self.assertEqual(abn.breakpoints[1].arm, "q")
        self.assertEqual(abn.breakpoints[1].region, 4)
        self.assertEqual(abn.breakpoints[1].band, 2)

    def test_parse_ring_x_chromosome(self):
        """Test r(X) - ring X chromosome."""
        result = self.parser.parse("45,X,r(X)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "r")
        self.assertEqual(abn.chromosome, "X")


class TestParserMarkerChromosomes(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_marker_single(self):
        """Test +mar - single marker chromosome."""
        result = self.parser.parse("47,XX,+mar")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "+mar")
        self.assertEqual(abn.chromosome, "mar")

    def test_parse_marker_multiple(self):
        """Test +2mar - two marker chromosomes."""
        result = self.parser.parse("48,XY,+2mar")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "+mar")
        self.assertEqual(abn.copy_count, 2)

    def test_parse_marker_numbered(self):
        """Test +mar1 - numbered marker."""
        result = self.parser.parse("47,XX,+mar1")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "+mar")
        self.assertEqual(abn.chromosome, "mar1")


class TestParserDerivativeChromosomes(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_derivative_translocation(self):
        """Test der(22)t(9;22)(q34;q11.2) - Philadelphia chromosome."""
        result = self.parser.parse("46,XX,der(22)t(9;22)(q34;q11.2)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "der")
        self.assertEqual(abn.chromosome, "22")
        self.assertIn("t(9;22)", abn.raw)

    def test_parse_derivative_deletion(self):
        """Test der(1)del(1)(p31) - derivative from deletion."""
        result = self.parser.parse("46,XY,der(1)del(1)(p31)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "der")
        self.assertEqual(abn.chromosome, "1")
        self.assertIn("del(1)", abn.raw)

    def test_parse_derivative_sex_chromosome(self):
        """Test der(X) - derivative X chromosome."""
        result = self.parser.parse("46,X,der(X)t(X;8)(p22;q24)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "der")
        self.assertEqual(abn.chromosome, "X")


class TestParserDoubleMinutesAndHSR(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_dmin(self):
        """Test dmin - double minutes."""
        result = self.parser.parse("47,XX,+21,dmin")
        # Find dmin abnormality
        dmin_abn = [a for a in result.abnormalities if a.type == "dmin"][0]
        self.assertEqual(dmin_abn.type, "dmin")

    def test_parse_hsr(self):
        """Test hsr - homogeneously staining region."""
        result = self.parser.parse("46,XX,hsr")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "hsr")

    def test_parse_hsr_with_location(self):
        """Test hsr(1)(p22) - HSR at specific location."""
        result = self.parser.parse("46,XY,hsr(1)(p22)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "hsr")
        self.assertEqual(abn.chromosome, "1")
        self.assertEqual(len(abn.breakpoints), 1)
        self.assertEqual(abn.breakpoints[0].arm, "p")


class TestParserMosaicism(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_mosaic_two_cell_lines(self):
        """Test 47,XX,+21[10]/46,XX[20] - mosaic with counts."""
        result = self.parser.parse("47,XX,+21[10]/46,XX[20]")
        self.assertIsNotNone(result.cell_lines)
        self.assertEqual(len(result.cell_lines), 2)
        # First cell line: 47,XX,+21[10]
        self.assertEqual(result.cell_lines[0].chromosome_count, 47)
        self.assertEqual(result.cell_lines[0].sex_chromosomes, "XX")
        self.assertEqual(len(result.cell_lines[0].abnormalities), 1)
        self.assertEqual(result.cell_lines[0].count, 10)
        # Second cell line: 46,XX[20]
        self.assertEqual(result.cell_lines[1].chromosome_count, 46)
        self.assertEqual(result.cell_lines[1].sex_chromosomes, "XX")
        self.assertEqual(result.cell_lines[1].abnormalities, [])
        self.assertEqual(result.cell_lines[1].count, 20)

    def test_parse_mosaic_without_counts(self):
        """Test 47,XX,+21/46,XX - mosaic without cell counts."""
        result = self.parser.parse("47,XX,+21/46,XX")
        self.assertIsNotNone(result.cell_lines)
        self.assertEqual(len(result.cell_lines), 2)
        self.assertEqual(result.cell_lines[0].count, 0)  # No count specified
        self.assertEqual(result.cell_lines[1].count, 0)

    def test_parse_mosaic_three_lines(self):
        """Test mosaic with three cell lines."""
        result = self.parser.parse("47,XX,+21[5]/46,XX,del(5)(q13)[10]/46,XX[15]")
        self.assertIsNotNone(result.cell_lines)
        self.assertEqual(len(result.cell_lines), 3)
        self.assertEqual(result.cell_lines[0].count, 5)
        self.assertEqual(result.cell_lines[1].count, 10)
        self.assertEqual(result.cell_lines[2].count, 15)
        # Second line has deletion
        self.assertEqual(len(result.cell_lines[1].abnormalities), 1)
        self.assertEqual(result.cell_lines[1].abnormalities[0].type, "del")


class TestParserUncertainty(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_uncertain_numerical(self):
        """Test ?+21 - uncertain trisomy."""
        result = self.parser.parse("47,XX,?+21")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "+")
        self.assertEqual(abn.chromosome, "21")
        self.assertTrue(abn.uncertain)

    def test_parse_uncertain_deletion(self):
        """Test ?del(5)(q13) - uncertain deletion."""
        result = self.parser.parse("46,XX,?del(5)(q13)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "del")
        self.assertTrue(abn.uncertain)

    def test_parse_uncertain_translocation(self):
        """Test ?t(9;22)(q34;q11) - uncertain translocation."""
        result = self.parser.parse("46,XX,?t(9;22)(q34;q11)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "t")
        self.assertTrue(abn.uncertain)

    def test_parse_certain_abnormality(self):
        """Test del(5)(q13) - certain deletion (no ?)."""
        result = self.parser.parse("46,XX,del(5)(q13)")
        abn = result.abnormalities[0]
        self.assertFalse(abn.uncertain)


class TestParserInheritance(unittest.TestCase):
    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_maternal_deletion(self):
        """Test del(5)(q13)mat - maternal inheritance."""
        result = self.parser.parse("46,XX,del(5)(q13)mat")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "del")
        self.assertEqual(abn.inheritance, "mat")

    def test_parse_paternal_translocation(self):
        """Test t(9;22)(q34;q11)pat - paternal inheritance."""
        result = self.parser.parse("46,XY,t(9;22)(q34;q11)pat")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "t")
        self.assertEqual(abn.inheritance, "pat")

    def test_parse_de_novo_trisomy(self):
        """Test +21dn - de novo trisomy."""
        result = self.parser.parse("47,XX,+21dn")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "+")
        self.assertEqual(abn.chromosome, "21")
        self.assertEqual(abn.inheritance, "dn")

    def test_parse_no_inheritance(self):
        """Test del(5)(q13) - no inheritance specified."""
        result = self.parser.parse("46,XX,del(5)(q13)")
        abn = result.abnormalities[0]
        self.assertIsNone(abn.inheritance)

    def test_parse_uncertain_with_inheritance(self):
        """Test ?del(5)(q13)mat - uncertain maternal."""
        result = self.parser.parse("46,XX,?del(5)(q13)mat")
        abn = result.abnormalities[0]
        self.assertTrue(abn.uncertain)
        self.assertEqual(abn.inheritance, "mat")


class TestParserInsertion(unittest.TestCase):
    """Tests for insertion abnormality parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_interchromosomal_insertion(self):
        """Test ins(5;2)(p14;q21q31) - segment from chr 2 inserted into chr 5."""
        result = self.parser.parse("46,XX,ins(5;2)(p14;q21q31)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "ins")
        self.assertEqual(abn.chromosome, "5;2")
        # Should have 3 breakpoints: insertion site and two segment ends
        self.assertEqual(len(abn.breakpoints), 3)

    def test_parse_intrachromosomal_insertion(self):
        """Test ins(2)(p13q21q31) - direct insertion within same chromosome."""
        result = self.parser.parse("46,XY,ins(2)(p13q21q31)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "ins")
        self.assertEqual(abn.chromosome, "2")
        # Should have 3 breakpoints
        self.assertEqual(len(abn.breakpoints), 3)

    def test_parse_insertion_x_chromosome(self):
        """Test insertion involving X chromosome."""
        result = self.parser.parse("46,XX,ins(X;1)(q21;p31p36)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "ins")
        self.assertEqual(abn.chromosome, "X;1")

    def test_parse_insertion_with_inheritance(self):
        """Test ins(5;2)(p14;q21q31)mat - maternal insertion."""
        result = self.parser.parse("46,XX,ins(5;2)(p14;q21q31)mat")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "ins")
        self.assertEqual(abn.inheritance, "mat")


class TestParserAdditionalMaterial(unittest.TestCase):
    """Tests for additional material (add) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_add_autosome(self):
        """Test add(7)(p22) - additional material on chr 7."""
        result = self.parser.parse("46,XX,add(7)(p22)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "add")
        self.assertEqual(abn.chromosome, "7")
        self.assertEqual(len(abn.breakpoints), 1)
        self.assertEqual(abn.breakpoints[0].arm, "p")

    def test_parse_add_x_chromosome(self):
        """Test add(X)(q28) - additional material on X."""
        result = self.parser.parse("46,XX,add(X)(q28)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "add")
        self.assertEqual(abn.chromosome, "X")

    def test_parse_add_with_uncertainty(self):
        """Test ?add(7)(p22) - uncertain additional material."""
        result = self.parser.parse("46,XY,?add(7)(p22)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "add")
        self.assertTrue(abn.uncertain)

    def test_parse_add_with_inheritance(self):
        """Test add(7)(p22)dn - de novo additional material."""
        result = self.parser.parse("46,XY,add(7)(p22)dn")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "add")
        self.assertEqual(abn.inheritance, "dn")


class TestParserTriplication(unittest.TestCase):
    """Tests for triplication (trp) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_triplication(self):
        """Test trp(1)(q21q32) - triplication of segment."""
        result = self.parser.parse("46,XX,trp(1)(q21q32)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "trp")
        self.assertEqual(abn.chromosome, "1")
        self.assertEqual(len(abn.breakpoints), 2)

    def test_parse_triplication_x_chromosome(self):
        """Test trp(X)(q21q28) - triplication on X."""
        result = self.parser.parse("46,XX,trp(X)(q21q28)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "trp")
        self.assertEqual(abn.chromosome, "X")

    def test_parse_triplication_with_inheritance(self):
        """Test trp(1)(q21q32)mat - maternal triplication."""
        result = self.parser.parse("46,XY,trp(1)(q21q32)mat")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "trp")
        self.assertEqual(abn.inheritance, "mat")


class TestParserDicentric(unittest.TestCase):
    """Tests for dicentric chromosome (dic) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_dicentric_two_chromosomes(self):
        """Test dic(13;14)(q14;q11) - dicentric from chr 13 and 14."""
        result = self.parser.parse("45,XX,dic(13;14)(q14;q11)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "dic")
        self.assertEqual(abn.chromosome, "13;14")
        self.assertEqual(len(abn.breakpoints), 2)

    def test_parse_dicentric_x_chromosome(self):
        """Test dic(X;Y)(p22;q11) - dicentric involving sex chromosomes."""
        result = self.parser.parse("45,X,dic(X;Y)(p22;q11)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "dic")
        self.assertEqual(abn.chromosome, "X;Y")

    def test_parse_dicentric_with_inheritance(self):
        """Test dic(13;14)(q14;q11)mat - maternal dicentric."""
        result = self.parser.parse("45,XX,dic(13;14)(q14;q11)mat")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "dic")
        self.assertEqual(abn.inheritance, "mat")


class TestParserIsodicentric(unittest.TestCase):
    """Tests for isodicentric chromosome (idic) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_isodicentric_y(self):
        """Test idic(Y)(q11) - isodicentric Y chromosome."""
        result = self.parser.parse("46,X,idic(Y)(q11)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "idic")
        self.assertEqual(abn.chromosome, "Y")
        self.assertEqual(len(abn.breakpoints), 1)

    def test_parse_isodicentric_autosome(self):
        """Test idic(15)(q13) - isodicentric chr 15."""
        result = self.parser.parse("46,XX,idic(15)(q13)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "idic")
        self.assertEqual(abn.chromosome, "15")

    def test_parse_isodicentric_with_inheritance(self):
        """Test idic(Y)(q11)dn - de novo isodicentric Y."""
        result = self.parser.parse("46,X,idic(Y)(q11)dn")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "idic")
        self.assertEqual(abn.inheritance, "dn")


class TestParserFragileSite(unittest.TestCase):
    """Tests for fragile site (fra) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_fragile_site_x(self):
        """Test fra(X)(q27.3) - fragile X site."""
        result = self.parser.parse("46,Y,fra(X)(q27.3)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "fra")
        self.assertEqual(abn.chromosome, "X")
        self.assertEqual(len(abn.breakpoints), 1)
        self.assertEqual(abn.breakpoints[0].arm, "q")

    def test_parse_fragile_site_autosome(self):
        """Test fra(16)(p13.11) - autosomal fragile site."""
        result = self.parser.parse("46,XX,fra(16)(p13.11)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "fra")
        self.assertEqual(abn.chromosome, "16")

    def test_parse_fragile_site_with_inheritance(self):
        """Test fra(X)(q27.3)mat - maternal fragile X."""
        result = self.parser.parse("46,Y,fra(X)(q27.3)mat")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "fra")
        self.assertEqual(abn.inheritance, "mat")


class TestParserRobertsonianTranslocation(unittest.TestCase):
    """Tests for Robertsonian translocation (rob) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_rob_13_14(self):
        """Test rob(13;14)(q10;q10) - common Robertsonian translocation."""
        result = self.parser.parse("45,XX,rob(13;14)(q10;q10)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "rob")
        self.assertEqual(abn.chromosome, "13;14")
        self.assertEqual(len(abn.breakpoints), 2)

    def test_parse_rob_14_21(self):
        """Test rob(14;21)(q10;q10) - another common Robertsonian."""
        result = self.parser.parse("45,XY,rob(14;21)(q10;q10)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "rob")
        self.assertEqual(abn.chromosome, "14;21")

    def test_parse_rob_with_inheritance(self):
        """Test rob(13;14)(q10;q10)mat - maternal Robertsonian."""
        result = self.parser.parse("45,XX,rob(13;14)(q10;q10)mat")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "rob")
        self.assertEqual(abn.inheritance, "mat")


class TestParserQuadruplication(unittest.TestCase):
    """Tests for quadruplication (qdp) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_quadruplication(self):
        """Test qdp(1)(q21q32) - quadruplication of segment."""
        result = self.parser.parse("46,XX,qdp(1)(q21q32)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "qdp")
        self.assertEqual(abn.chromosome, "1")
        self.assertEqual(len(abn.breakpoints), 2)

    def test_parse_quadruplication_x_chromosome(self):
        """Test qdp(X)(p21p11) - quadruplication on X."""
        result = self.parser.parse("46,XX,qdp(X)(p21p11)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "qdp")
        self.assertEqual(abn.chromosome, "X")

    def test_parse_quadruplication_with_inheritance(self):
        """Test qdp(1)(q21q32)dn - de novo quadruplication."""
        result = self.parser.parse("46,XY,qdp(1)(q21q32)dn")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "qdp")
        self.assertEqual(abn.inheritance, "dn")


class TestParserPseudodicentric(unittest.TestCase):
    """Tests for pseudodicentric chromosome (psu dic) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_pseudodicentric_two_chromosomes(self):
        """Test psu dic(13;14)(q14;q11) - pseudodicentric from chr 13 and 14."""
        result = self.parser.parse("45,XX,psu dic(13;14)(q14;q11)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "psu dic")
        self.assertEqual(abn.chromosome, "13;14")
        self.assertEqual(len(abn.breakpoints), 2)

    def test_parse_pseudodicentric_no_space(self):
        """Test psudic(13;14)(q14;q11) - pseudodicentric without space."""
        result = self.parser.parse("45,XX,psu dic(13;14)(q14;q11)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "psu dic")
        self.assertEqual(abn.chromosome, "13;14")

    def test_parse_pseudodicentric_with_inheritance(self):
        """Test psu dic(13;14)(q14;q11)mat - maternal pseudodicentric."""
        result = self.parser.parse("45,XX,psu dic(13;14)(q14;q11)mat")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "psu dic")
        self.assertEqual(abn.inheritance, "mat")


class TestParserAcentricFragment(unittest.TestCase):
    """Tests for acentric fragment (ace) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_acentric_fragment(self):
        """Test ace(1)(q21q31) - acentric fragment from chr 1."""
        result = self.parser.parse("46,XX,ace(1)(q21q31)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "ace")
        self.assertEqual(abn.chromosome, "1")
        self.assertEqual(len(abn.breakpoints), 2)

    def test_parse_acentric_single_breakpoint(self):
        """Test ace(5)(p15) - acentric with single breakpoint."""
        result = self.parser.parse("46,XY,ace(5)(p15)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "ace")
        self.assertEqual(abn.chromosome, "5")
        self.assertEqual(len(abn.breakpoints), 1)

    def test_parse_acentric_x_chromosome(self):
        """Test ace(X)(q21q28) - acentric X chromosome fragment."""
        result = self.parser.parse("46,XX,ace(X)(q21q28)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "ace")
        self.assertEqual(abn.chromosome, "X")

    def test_parse_acentric_with_uncertainty(self):
        """Test ?ace(1)(q21q31) - uncertain acentric."""
        result = self.parser.parse("46,XX,?ace(1)(q21q31)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "ace")
        self.assertTrue(abn.uncertain)


class TestParserTelomericAssociation(unittest.TestCase):
    """Tests for telomeric association (tas) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_telomeric_association(self):
        """Test tas(13;14)(p11;p11) - telomeric association."""
        result = self.parser.parse("46,XX,tas(13;14)(p11;p11)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "tas")
        self.assertEqual(abn.chromosome, "13;14")
        self.assertEqual(len(abn.breakpoints), 2)

    def test_parse_telomeric_association_three_chromosomes(self):
        """Test tas(13;14;15)(p11;p11;p11) - three-way telomeric association."""
        result = self.parser.parse("46,XY,tas(13;14;15)(p11;p11;p11)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "tas")
        self.assertEqual(abn.chromosome, "13;14;15")
        self.assertEqual(len(abn.breakpoints), 3)

    def test_parse_telomeric_association_with_inheritance(self):
        """Test tas(13;14)(p11;p11)dn - de novo telomeric association."""
        result = self.parser.parse("46,XX,tas(13;14)(p11;p11)dn")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "tas")
        self.assertEqual(abn.inheritance, "dn")


class TestParserFission(unittest.TestCase):
    """Tests for fission (fis) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_fission_p_arm(self):
        """Test fis(1)(p10) - fission at p arm centromere."""
        result = self.parser.parse("47,XX,fis(1)(p10)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "fis")
        self.assertEqual(abn.chromosome, "1")
        self.assertEqual(len(abn.breakpoints), 1)
        self.assertEqual(abn.breakpoints[0].arm, "p")

    def test_parse_fission_q_arm(self):
        """Test fis(1)(q10) - fission at q arm centromere."""
        result = self.parser.parse("47,XY,fis(1)(q10)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "fis")
        self.assertEqual(abn.breakpoints[0].arm, "q")

    def test_parse_fission_x_chromosome(self):
        """Test fis(X)(q10) - fission of X chromosome."""
        result = self.parser.parse("47,XX,fis(X)(q10)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "fis")
        self.assertEqual(abn.chromosome, "X")

    def test_parse_fission_with_inheritance(self):
        """Test fis(1)(p10)dn - de novo fission."""
        result = self.parser.parse("47,XX,fis(1)(p10)dn")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "fis")
        self.assertEqual(abn.inheritance, "dn")


class TestParserNeocentromere(unittest.TestCase):
    """Tests for neocentromere (neo) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_neocentromere(self):
        """Test neo(1)(q21) - neocentromere at chr 1 q21."""
        result = self.parser.parse("46,XX,neo(1)(q21)")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "neo")
        self.assertEqual(abn.chromosome, "1")
        self.assertEqual(len(abn.breakpoints), 1)
        self.assertEqual(abn.breakpoints[0].arm, "q")

    def test_parse_neocentromere_x_chromosome(self):
        """Test neo(X)(p11) - neocentromere on X chromosome."""
        result = self.parser.parse("46,XX,neo(X)(p11)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "neo")
        self.assertEqual(abn.chromosome, "X")

    def test_parse_neocentromere_with_subband(self):
        """Test neo(15)(q24.2) - neocentromere with subband."""
        result = self.parser.parse("46,XY,neo(15)(q24.2)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "neo")
        self.assertEqual(abn.breakpoints[0].subband, "2")

    def test_parse_neocentromere_with_uncertainty(self):
        """Test ?neo(1)(q21) - uncertain neocentromere."""
        result = self.parser.parse("46,XX,?neo(1)(q21)")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "neo")
        self.assertTrue(abn.uncertain)


class TestParserIncomplete(unittest.TestCase):
    """Tests for incomplete karyotype (inc) parsing."""

    def setUp(self):
        self.parser = KaryotypeParser()

    def test_parse_incomplete_simple(self):
        """Test 46,XX,inc - simple incomplete karyotype."""
        result = self.parser.parse("46,XX,inc")
        self.assertEqual(len(result.abnormalities), 1)
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "inc")
        self.assertEqual(abn.chromosome, "")
        self.assertEqual(len(abn.breakpoints), 0)

    def test_parse_incomplete_with_other_abnormality(self):
        """Test 47,XX,+21,inc - incomplete with other abnormality."""
        result = self.parser.parse("47,XX,+21,inc")
        self.assertEqual(len(result.abnormalities), 2)
        inc_abn = [a for a in result.abnormalities if a.type == "inc"][0]
        self.assertEqual(inc_abn.type, "inc")

    def test_parse_incomplete_male(self):
        """Test 46,XY,inc - incomplete male karyotype."""
        result = self.parser.parse("46,XY,inc")
        abn = result.abnormalities[0]
        self.assertEqual(abn.type, "inc")


if __name__ == "__main__":
    unittest.main()
