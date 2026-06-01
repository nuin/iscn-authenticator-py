# iscn_authenticator/parser.py
"""Parser for ISCN karyotype strings."""

import re

from iscn_authenticator.models import Abnormality, Breakpoint, CellLine, KaryotypeAST

# Pattern for cell line count: [10], [20], etc.
CELL_LINE_COUNT_PATTERN = re.compile(r"^(.+?)\[(\d+)\]$")


class ParseError(Exception):
    """Raised when karyotype string cannot be parsed."""

    pass


class KaryotypeParser:
    """Parses ISCN karyotype strings into AST."""

    # Regex patterns
    SEX_CHROMOSOMES_PATTERN = re.compile(r"^[XYU]+$")
    NUMERICAL_ABNORMALITY_PATTERN = re.compile(r"^([+-])(\d{1,2}|[XY])$")
    DELETION_PATTERN = re.compile(r"^del\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    DUPLICATION_PATTERN = re.compile(r"^dup\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    INVERSION_PATTERN = re.compile(r"^inv\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    TRANSLOCATION_PATTERN = re.compile(r"^t\(([^)]+)\)\(([^)]+)\)$")
    BREAKPOINT_PATTERN = re.compile(r"^([pq])(\d+)(?:\.(\d+))?$")
    # Isochromosome: i(17q) short form or i(17)(q10) long form
    ISOCHROMOSOME_SHORT_PATTERN = re.compile(r"^i\((\d{1,2}|[XY])([pq])\)$")
    ISOCHROMOSOME_LONG_PATTERN = re.compile(r"^i\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    # Ring chromosome: r(1) simple or r(1)(p36q42) with breakpoints
    RING_SIMPLE_PATTERN = re.compile(r"^r\((\d{1,2}|[XY])\)$")
    RING_BREAKPOINT_PATTERN = re.compile(r"^r\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    # Marker chromosome: +mar, +2mar, +mar1
    MARKER_PATTERN = re.compile(r"^\+(\d*)mar(\d*)$")
    # Derivative chromosome: der(22)t(9;22)(...) or der(1)del(1)(...)
    DERIVATIVE_PATTERN = re.compile(r"^der\((\d{1,2}|[XY])\)(.+)$")
    # Double minutes: dmin
    DMIN_PATTERN = re.compile(r"^dmin$")
    # Homogeneously staining region: hsr or hsr(1)(p22)
    HSR_SIMPLE_PATTERN = re.compile(r"^hsr$")
    HSR_LOCATION_PATTERN = re.compile(r"^hsr\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    # Insertion: ins(5;2)(p14;q21q31) or ins(2)(p13q21q31)
    INSERTION_PATTERN = re.compile(r"^ins\(([^)]+)\)\(([^)]+)\)$")
    # Additional material of unknown origin: add(7)(p22)
    ADD_PATTERN = re.compile(r"^add\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    # Triplication: trp(1)(q21q32)
    TRIPLICATION_PATTERN = re.compile(r"^trp\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    # Dicentric chromosome: dic(13;14)(q14;q11)
    DICENTRIC_PATTERN = re.compile(r"^dic\(([^)]+)\)\(([^)]+)\)$")
    # Isodicentric chromosome: idic(Y)(q11)
    ISODICENTRIC_PATTERN = re.compile(r"^idic\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    # Fragile site: fra(X)(q27.3)
    FRAGILE_SITE_PATTERN = re.compile(r"^fra\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    # Robertsonian translocation: rob(13;14)(q10;q10)
    ROBERTSONIAN_PATTERN = re.compile(r"^rob\(([^)]+)\)\(([^)]+)\)$")
    # Quadruplication: qdp(1)(q21q32)
    QUADRUPLICATION_PATTERN = re.compile(r"^qdp\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    # Pseudodicentric: psu dic(13;14)(q14;q11)
    PSEUDODICENTRIC_PATTERN = re.compile(r"^psu\s*dic\(([^)]+)\)\(([^)]+)\)$")
    # Acentric fragment: ace(1)(q21q31)
    ACENTRIC_PATTERN = re.compile(r"^ace\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    # Telomeric association: tas(13;14)(p11;p11)
    TELOMERIC_ASSOC_PATTERN = re.compile(r"^tas\(([^)]+)\)\(([^)]+)\)$")
    # Fission: fis(1)(p10) or fis(1)(q10)
    FISSION_PATTERN = re.compile(r"^fis\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    # Neocentromere: neo(1)(q21)
    NEOCENTROMERE_PATTERN = re.compile(r"^neo\((\d{1,2}|[XY])\)\(([^)]+)\)$")
    # Incomplete karyotype: inc
    INCOMPLETE_PATTERN = re.compile(r"^inc$")

    def parse(self, karyotype: str) -> KaryotypeAST:
        """Parse a karyotype string into an AST."""
        if not karyotype or not karyotype.strip():
            raise ParseError("Karyotype string is empty")

        karyotype = karyotype.strip()

        # Check for mosaicism (cell lines separated by /)
        if "/" in karyotype:
            return self._parse_mosaic(karyotype)

        return self._parse_single_karyotype(karyotype)

    def _parse_single_karyotype(
        self, karyotype: str, extract_count: bool = False
    ) -> tuple[KaryotypeAST, int] | KaryotypeAST:
        """Parse a single karyotype string (non-mosaic)."""
        count = 0

        # Extract cell count if present (e.g., "46,XX[10]")
        if extract_count:
            count_match = CELL_LINE_COUNT_PATTERN.match(karyotype)
            if count_match:
                karyotype = count_match.group(1)
                count = int(count_match.group(2))

        # Split on comma
        if "," not in karyotype:
            raise ParseError("Missing comma separator between chromosome count and sex chromosomes")

        parts = karyotype.split(",")

        chromosome_count = self._parse_chromosome_count(parts[0])

        sex_chromosomes = self._parse_sex_chromosomes(parts[1])

        # Parse abnormalities (if any)
        abnormalities = []
        if len(parts) > 2:
            abnormalities = self._parse_abnormalities(parts[2:])

        ast = KaryotypeAST(
            chromosome_count=chromosome_count,
            sex_chromosomes=sex_chromosomes,
            abnormalities=abnormalities,
            cell_lines=None,
            modifiers=None,
        )

        if extract_count:
            return ast, count
        return ast

    def _parse_mosaic(self, karyotype: str) -> KaryotypeAST:
        """Parse a mosaic karyotype with multiple cell lines."""
        cell_line_strs = karyotype.split("/")
        cell_lines = []

        for line_str in cell_line_strs:
            line_str = line_str.strip()
            ast, count = self._parse_single_karyotype(line_str, extract_count=True)
            cell_line = CellLine(
                chromosome_count=ast.chromosome_count,
                sex_chromosomes=ast.sex_chromosomes,
                abnormalities=ast.abnormalities,
                count=count,
                is_donor=False,
            )
            cell_lines.append(cell_line)

        # Use the first cell line as the main karyotype info
        first = cell_lines[0]
        return KaryotypeAST(
            chromosome_count=first.chromosome_count,
            sex_chromosomes=first.sex_chromosomes,
            abnormalities=first.abnormalities,
            cell_lines=cell_lines,
            modifiers=None,
        )

    def _parse_chromosome_count(self, count_str: str) -> int | str:
        """Parse chromosome count (number or range)."""
        count_str = count_str.strip()

        # Handle range notation (e.g., "45~48")
        if "~" in count_str:
            return count_str

        if not count_str.isdigit():
            raise ParseError(f"Invalid chromosome count: '{count_str}' is not a number")

        return int(count_str)

    def _parse_sex_chromosomes(self, sex_str: str) -> str:
        """Parse sex chromosome designation."""
        sex_str = sex_str.strip()

        if not self.SEX_CHROMOSOMES_PATTERN.match(sex_str):
            raise ParseError(f"Invalid sex chromosomes: '{sex_str}' must contain only X, Y, or U")

        return sex_str

    def _parse_breakpoint(self, bp_str: str) -> Breakpoint:
        """Parse a single breakpoint like 'q13' or 'p11.2'."""
        match = self.BREAKPOINT_PATTERN.match(bp_str)
        if not match:
            raise ParseError(f"Invalid breakpoint format: '{bp_str}'")

        arm = match.group(1)
        region_band = match.group(2)
        subband = match.group(3)

        # Split region and band (e.g., "13" -> region=1, band=3)
        if len(region_band) >= 2:
            region = int(region_band[0])
            band = int(region_band[1:])
        else:
            region = int(region_band)
            band = 0

        return Breakpoint(arm=arm, region=region, band=band, subband=subband, uncertain=False)

    def _parse_deletion(self, part: str) -> Abnormality:
        """Parse a deletion abnormality."""
        match = self.DELETION_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid deletion format: '{part}'")

        chromosome = match.group(1)
        breakpoint_str = match.group(2)

        # Parse breakpoints (could be single or double)
        breakpoints = []
        # Check for interstitial deletion (two breakpoints like q13q33)
        double_bp = re.match(r"^([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)$", breakpoint_str)
        if double_bp:
            breakpoints.append(self._parse_breakpoint(double_bp.group(1)))
            breakpoints.append(self._parse_breakpoint(double_bp.group(2)))
        else:
            # Single breakpoint (terminal deletion)
            breakpoints.append(self._parse_breakpoint(breakpoint_str))

        return Abnormality(
            type="del",
            chromosome=chromosome,
            breakpoints=breakpoints,
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_duplication(self, part: str) -> Abnormality:
        """Parse a duplication abnormality."""
        match = self.DUPLICATION_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid duplication format: '{part}'")

        chromosome = match.group(1)
        breakpoint_str = match.group(2)

        # Parse breakpoints (could be single or double)
        breakpoints = []
        double_bp = re.match(r"^([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)$", breakpoint_str)
        if double_bp:
            breakpoints.append(self._parse_breakpoint(double_bp.group(1)))
            breakpoints.append(self._parse_breakpoint(double_bp.group(2)))
        else:
            breakpoints.append(self._parse_breakpoint(breakpoint_str))

        return Abnormality(
            type="dup",
            chromosome=chromosome,
            breakpoints=breakpoints,
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_inversion(self, part: str) -> Abnormality:
        """Parse an inversion abnormality."""
        match = self.INVERSION_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid inversion format: '{part}'")

        chromosome = match.group(1)
        breakpoint_str = match.group(2)

        # Inversions always have two breakpoints
        breakpoints = []
        double_bp = re.match(r"^([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)$", breakpoint_str)
        if double_bp:
            breakpoints.append(self._parse_breakpoint(double_bp.group(1)))
            breakpoints.append(self._parse_breakpoint(double_bp.group(2)))
        else:
            raise ParseError(f"Inversion requires two breakpoints: '{part}'")

        return Abnormality(
            type="inv",
            chromosome=chromosome,
            breakpoints=breakpoints,
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_translocation(self, part: str) -> Abnormality:
        """Parse a translocation abnormality."""
        match = self.TRANSLOCATION_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid translocation format: '{part}'")

        chromosomes_str = match.group(1)  # e.g., "9;22" or "1;3;5"
        breakpoints_str = match.group(2)  # e.g., "q34;q11.2" or "p32;q21;q31"

        # Parse breakpoints (semicolon-separated)
        bp_parts = breakpoints_str.split(";")
        breakpoints = [self._parse_breakpoint(bp.strip()) for bp in bp_parts]

        return Abnormality(
            type="t",
            chromosome=chromosomes_str,
            breakpoints=breakpoints,
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_isochromosome(self, part: str) -> Abnormality:
        """Parse an isochromosome abnormality."""
        # Try short form first: i(17q) or i(Xq)
        short_match = self.ISOCHROMOSOME_SHORT_PATTERN.match(part)
        if short_match:
            chromosome = short_match.group(1)
            arm = short_match.group(2)
            # Create a breakpoint with just the arm (region/band at centromere)
            breakpoint = Breakpoint(arm=arm, region=1, band=0, subband=None, uncertain=False)
            return Abnormality(
                type="i",
                chromosome=chromosome,
                breakpoints=[breakpoint],
                inheritance=None,
                uncertain=False,
                copy_count=None,
                raw=part,
            )

        # Try long form: i(17)(q10)
        long_match = self.ISOCHROMOSOME_LONG_PATTERN.match(part)
        if long_match:
            chromosome = long_match.group(1)
            breakpoint_str = long_match.group(2)
            breakpoint = self._parse_breakpoint(breakpoint_str)
            return Abnormality(
                type="i",
                chromosome=chromosome,
                breakpoints=[breakpoint],
                inheritance=None,
                uncertain=False,
                copy_count=None,
                raw=part,
            )

        raise ParseError(f"Invalid isochromosome format: '{part}'")

    def _parse_ring(self, part: str) -> Abnormality:
        """Parse a ring chromosome abnormality."""
        # Try simple form first: r(1)
        simple_match = self.RING_SIMPLE_PATTERN.match(part)
        if simple_match:
            chromosome = simple_match.group(1)
            return Abnormality(
                type="r",
                chromosome=chromosome,
                breakpoints=[],
                inheritance=None,
                uncertain=False,
                copy_count=None,
                raw=part,
            )

        # Try breakpoint form: r(1)(p36q42)
        bp_match = self.RING_BREAKPOINT_PATTERN.match(part)
        if bp_match:
            chromosome = bp_match.group(1)
            breakpoint_str = bp_match.group(2)
            # Parse two breakpoints (p arm and q arm)
            breakpoints = []
            double_bp = re.match(r"^([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)$", breakpoint_str)
            if double_bp:
                breakpoints.append(self._parse_breakpoint(double_bp.group(1)))
                breakpoints.append(self._parse_breakpoint(double_bp.group(2)))
            return Abnormality(
                type="r",
                chromosome=chromosome,
                breakpoints=breakpoints,
                inheritance=None,
                uncertain=False,
                copy_count=None,
                raw=part,
            )

        raise ParseError(f"Invalid ring chromosome format: '{part}'")

    def _parse_insertion(self, part: str) -> Abnormality:
        """Parse an insertion abnormality.

        Formats:
        - ins(5;2)(p14;q21q31) - interchromosomal: segment from chr 2 inserted into chr 5
        - ins(2)(p13q21q31) - intrachromosomal: direct insertion within same chromosome
        """
        match = self.INSERTION_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid insertion format: '{part}'")

        chromosomes_str = match.group(1)  # e.g., "5;2" or "2"
        breakpoints_str = match.group(2)  # e.g., "p14;q21q31" or "p13q21q31"

        breakpoints = []
        if ";" in breakpoints_str:
            # Interchromosomal: breakpoints separated by semicolon
            # Format: insertion_site;segment_start segment_end (e.g., "p14;q21q31")
            bp_parts = breakpoints_str.split(";")
            # First part is insertion site
            breakpoints.append(self._parse_breakpoint(bp_parts[0].strip()))
            # Second part contains two breakpoints (segment boundaries)
            segment_str = bp_parts[1].strip()
            double_bp = re.match(r"^([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)$", segment_str)
            if double_bp:
                breakpoints.append(self._parse_breakpoint(double_bp.group(1)))
                breakpoints.append(self._parse_breakpoint(double_bp.group(2)))
            else:
                # Single breakpoint after semicolon
                breakpoints.append(self._parse_breakpoint(segment_str))
        else:
            # Intrachromosomal: three consecutive breakpoints (e.g., "p13q21q31")
            # Try to parse three breakpoints
            triple_bp = re.match(
                r"^([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)$",
                breakpoints_str,
            )
            if triple_bp:
                breakpoints.append(self._parse_breakpoint(triple_bp.group(1)))
                breakpoints.append(self._parse_breakpoint(triple_bp.group(2)))
                breakpoints.append(self._parse_breakpoint(triple_bp.group(3)))
            else:
                raise ParseError(f"Invalid insertion breakpoints: '{breakpoints_str}'")

        return Abnormality(
            type="ins",
            chromosome=chromosomes_str,
            breakpoints=breakpoints,
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_add(self, part: str) -> Abnormality:
        """Parse additional material of unknown origin.

        Format: add(7)(p22) - additional material attached at chr 7 p22
        """
        match = self.ADD_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid additional material format: '{part}'")

        chromosome = match.group(1)
        breakpoint_str = match.group(2)
        breakpoint = self._parse_breakpoint(breakpoint_str)

        return Abnormality(
            type="add",
            chromosome=chromosome,
            breakpoints=[breakpoint],
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_triplication(self, part: str) -> Abnormality:
        """Parse a triplication abnormality.

        Format: trp(1)(q21q32) - triplication of segment q21 to q32 on chr 1
        """
        match = self.TRIPLICATION_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid triplication format: '{part}'")

        chromosome = match.group(1)
        breakpoint_str = match.group(2)

        # Parse two breakpoints (segment boundaries)
        breakpoints = []
        double_bp = re.match(r"^([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)$", breakpoint_str)
        if double_bp:
            breakpoints.append(self._parse_breakpoint(double_bp.group(1)))
            breakpoints.append(self._parse_breakpoint(double_bp.group(2)))
        else:
            raise ParseError(f"Triplication requires two breakpoints: '{part}'")

        return Abnormality(
            type="trp",
            chromosome=chromosome,
            breakpoints=breakpoints,
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_dicentric(self, part: str) -> Abnormality:
        """Parse a dicentric chromosome abnormality.

        Format: dic(13;14)(q14;q11) - dicentric formed from chr 13 and 14
        """
        match = self.DICENTRIC_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid dicentric format: '{part}'")

        chromosomes_str = match.group(1)  # e.g., "13;14"
        breakpoints_str = match.group(2)  # e.g., "q14;q11"

        # Parse breakpoints (semicolon-separated)
        bp_parts = breakpoints_str.split(";")
        breakpoints = [self._parse_breakpoint(bp.strip()) for bp in bp_parts]

        return Abnormality(
            type="dic",
            chromosome=chromosomes_str,
            breakpoints=breakpoints,
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_isodicentric(self, part: str) -> Abnormality:
        """Parse an isodicentric chromosome abnormality.

        Format: idic(Y)(q11) - isodicentric Y with breakpoint at q11
        """
        match = self.ISODICENTRIC_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid isodicentric format: '{part}'")

        chromosome = match.group(1)
        breakpoint_str = match.group(2)
        breakpoint = self._parse_breakpoint(breakpoint_str)

        return Abnormality(
            type="idic",
            chromosome=chromosome,
            breakpoints=[breakpoint],
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_fragile_site(self, part: str) -> Abnormality:
        """Parse a fragile site.

        Format: fra(X)(q27.3) - fragile site at Xq27.3
        """
        match = self.FRAGILE_SITE_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid fragile site format: '{part}'")

        chromosome = match.group(1)
        breakpoint_str = match.group(2)
        breakpoint = self._parse_breakpoint(breakpoint_str)

        return Abnormality(
            type="fra",
            chromosome=chromosome,
            breakpoints=[breakpoint],
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_robertsonian(self, part: str) -> Abnormality:
        """Parse a Robertsonian translocation.

        Format: rob(13;14)(q10;q10) - Robertsonian translocation between chr 13 and 14
        """
        match = self.ROBERTSONIAN_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid Robertsonian translocation format: '{part}'")

        chromosomes_str = match.group(1)  # e.g., "13;14"
        breakpoints_str = match.group(2)  # e.g., "q10;q10"

        # Parse breakpoints (semicolon-separated)
        bp_parts = breakpoints_str.split(";")
        breakpoints = [self._parse_breakpoint(bp.strip()) for bp in bp_parts]

        return Abnormality(
            type="rob",
            chromosome=chromosomes_str,
            breakpoints=breakpoints,
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_quadruplication(self, part: str) -> Abnormality:
        """Parse a quadruplication abnormality.

        Format: qdp(1)(q21q32) - quadruplication of segment q21 to q32 on chr 1
        """
        match = self.QUADRUPLICATION_PATTERN.match(part)
        if not match:
            raise ParseError(f"Invalid quadruplication format: '{part}'")

        chromosome = match.group(1)
        breakpoint_str = match.group(2)

        # Parse two breakpoints (segment boundaries)
        breakpoints = []
        double_bp = re.match(r"^([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)$", breakpoint_str)
        if double_bp:
            breakpoints.append(self._parse_breakpoint(double_bp.group(1)))
            breakpoints.append(self._parse_breakpoint(double_bp.group(2)))
        else:
            raise ParseError(f"Quadruplication requires two breakpoints: '{part}'")

        return Abnormality(
            type="qdp",
            chromosome=chromosome,
            breakpoints=breakpoints,
            inheritance=None,
            uncertain=False,
            copy_count=None,
            raw=part,
        )

    def _parse_abnormalities(self, parts: list[str]) -> list[Abnormality]:
        """Parse abnormality parts."""
        abnormalities = []
        for part in parts:
            part = part.strip()
            if not part:
                continue

            # Check for uncertainty marker (?)
            uncertain = False
            original_part = part
            if part.startswith("?"):
                uncertain = True
                part = part[1:]  # Remove the ? prefix

            # Check for inheritance notation (mat, pat, dn) at end
            inheritance = None
            if part.endswith("mat"):
                inheritance = "mat"
                part = part[:-3]
            elif part.endswith("pat"):
                inheritance = "pat"
                part = part[:-3]
            elif part.endswith("dn"):
                inheritance = "dn"
                part = part[:-2]

            # Try numerical abnormality (+21, -7, +X, -Y)
            num_match = self.NUMERICAL_ABNORMALITY_PATTERN.match(part)
            if num_match:
                abnormalities.append(
                    Abnormality(
                        type=num_match.group(1),  # "+" or "-"
                        chromosome=num_match.group(2),
                        breakpoints=[],
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=None,
                        raw=original_part,
                    )
                )
                continue

            # Try deletion
            if part.startswith("del("):
                abn = self._parse_deletion(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try additional material
            if part.startswith("add("):
                abn = self._parse_add(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try duplication
            if part.startswith("dup("):
                abn = self._parse_duplication(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try dicentric (must check before derivative since both start with 'd')
            if part.startswith("dic("):
                abn = self._parse_dicentric(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try isodicentric (must check before isochromosome since both start with 'i')
            if part.startswith("idic("):
                abn = self._parse_isodicentric(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try fragile site
            if part.startswith("fra("):
                abn = self._parse_fragile_site(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try inversion
            if part.startswith("inv("):
                abn = self._parse_inversion(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try triplication (must check before translocation since both start with 't')
            if part.startswith("trp("):
                abn = self._parse_triplication(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try quadruplication
            if part.startswith("qdp("):
                abn = self._parse_quadruplication(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try translocation
            if part.startswith("t("):
                abn = self._parse_translocation(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try insertion (must check before isochromosome since both start with 'i')
            if part.startswith("ins("):
                abn = self._parse_insertion(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try isochromosome
            if part.startswith("i("):
                abn = self._parse_isochromosome(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try Robertsonian translocation (must check before ring since both start with 'r')
            if part.startswith("rob("):
                abn = self._parse_robertsonian(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try ring chromosome
            if part.startswith("r("):
                abn = self._parse_ring(part)
                abn.uncertain = uncertain
                abn.inheritance = inheritance
                abn.raw = original_part
                abnormalities.append(abn)
                continue

            # Try marker chromosome (+mar, +2mar, +mar1)
            mar_match = self.MARKER_PATTERN.match(part)
            if mar_match:
                count_prefix = mar_match.group(1)  # e.g., "2" in +2mar
                marker_suffix = mar_match.group(2)  # e.g., "1" in +mar1
                copy_count = int(count_prefix) if count_prefix else None
                chromosome = "mar" + marker_suffix if marker_suffix else "mar"
                abnormalities.append(
                    Abnormality(
                        type="+mar",
                        chromosome=chromosome,
                        breakpoints=[],
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=copy_count,
                        raw=original_part,
                    )
                )
                continue

            # Try derivative chromosome (der(22)t(9;22)(...))
            der_match = self.DERIVATIVE_PATTERN.match(part)
            if der_match:
                chromosome = der_match.group(1)
                # The rearrangement description is captured in raw for now
                abnormalities.append(
                    Abnormality(
                        type="der",
                        chromosome=chromosome,
                        breakpoints=[],
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=None,
                        raw=original_part,
                    )
                )
                continue

            # Try double minutes (dmin)
            if self.DMIN_PATTERN.match(part):
                abnormalities.append(
                    Abnormality(
                        type="dmin",
                        chromosome="",
                        breakpoints=[],
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=None,
                        raw=original_part,
                    )
                )
                continue

            # Try HSR with location (hsr(1)(p22))
            hsr_loc_match = self.HSR_LOCATION_PATTERN.match(part)
            if hsr_loc_match:
                chromosome = hsr_loc_match.group(1)
                breakpoint_str = hsr_loc_match.group(2)
                breakpoint = self._parse_breakpoint(breakpoint_str)
                abnormalities.append(
                    Abnormality(
                        type="hsr",
                        chromosome=chromosome,
                        breakpoints=[breakpoint],
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=None,
                        raw=original_part,
                    )
                )
                continue

            # Try simple HSR (hsr)
            if self.HSR_SIMPLE_PATTERN.match(part):
                abnormalities.append(
                    Abnormality(
                        type="hsr",
                        chromosome="",
                        breakpoints=[],
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=None,
                        raw=original_part,
                    )
                )
                continue

            # Try pseudodicentric (psu dic(...))
            psu_match = self.PSEUDODICENTRIC_PATTERN.match(part)
            if psu_match:
                chromosomes = psu_match.group(1)
                breakpoint_str = psu_match.group(2)
                breakpoints = self._parse_multiple_breakpoints(breakpoint_str)
                abnormalities.append(
                    Abnormality(
                        type="psu dic",
                        chromosome=chromosomes,
                        breakpoints=breakpoints,
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=None,
                        raw=original_part,
                    )
                )
                continue

            # Try acentric fragment
            ace_match = self.ACENTRIC_PATTERN.match(part)
            if ace_match:
                chromosome = ace_match.group(1)
                breakpoint_str = ace_match.group(2)
                breakpoints = self._parse_breakpoints(breakpoint_str)
                abnormalities.append(
                    Abnormality(
                        type="ace",
                        chromosome=chromosome,
                        breakpoints=breakpoints,
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=None,
                        raw=original_part,
                    )
                )
                continue

            # Try telomeric association
            tas_match = self.TELOMERIC_ASSOC_PATTERN.match(part)
            if tas_match:
                chromosomes = tas_match.group(1)
                breakpoint_str = tas_match.group(2)
                breakpoints = self._parse_multiple_breakpoints(breakpoint_str)
                abnormalities.append(
                    Abnormality(
                        type="tas",
                        chromosome=chromosomes,
                        breakpoints=breakpoints,
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=None,
                        raw=original_part,
                    )
                )
                continue

            # Try fission
            fis_match = self.FISSION_PATTERN.match(part)
            if fis_match:
                chromosome = fis_match.group(1)
                breakpoint_str = fis_match.group(2)
                breakpoint = self._parse_breakpoint(breakpoint_str)
                abnormalities.append(
                    Abnormality(
                        type="fis",
                        chromosome=chromosome,
                        breakpoints=[breakpoint],
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=None,
                        raw=original_part,
                    )
                )
                continue

            # Try neocentromere
            neo_match = self.NEOCENTROMERE_PATTERN.match(part)
            if neo_match:
                chromosome = neo_match.group(1)
                breakpoint_str = neo_match.group(2)
                breakpoint = self._parse_breakpoint(breakpoint_str)
                abnormalities.append(
                    Abnormality(
                        type="neo",
                        chromosome=chromosome,
                        breakpoints=[breakpoint],
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=None,
                        raw=original_part,
                    )
                )
                continue

            # Try incomplete karyotype
            if self.INCOMPLETE_PATTERN.match(part):
                abnormalities.append(
                    Abnormality(
                        type="inc",
                        chromosome="",
                        breakpoints=[],
                        inheritance=inheritance,
                        uncertain=uncertain,
                        copy_count=None,
                        raw=original_part,
                    )
                )
                continue

            # Unknown abnormality type (will be expanded in later tasks)
            abnormalities.append(
                Abnormality(
                    type="unknown",
                    chromosome="",
                    breakpoints=[],
                    inheritance=inheritance,
                    uncertain=uncertain,
                    copy_count=None,
                    raw=original_part,
                )
            )
        return abnormalities

    def _parse_breakpoints(self, bp_str: str) -> list[Breakpoint]:
        """Parse one or two concatenated breakpoints like 'q21' or 'q21q31'."""
        breakpoints = []
        # Try to match two breakpoints
        double_bp = re.match(r"^([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)$", bp_str)
        if double_bp:
            breakpoints.append(self._parse_breakpoint(double_bp.group(1)))
            breakpoints.append(self._parse_breakpoint(double_bp.group(2)))
        else:
            # Single breakpoint
            breakpoints.append(self._parse_breakpoint(bp_str))
        return breakpoints

    def _parse_multiple_breakpoints(self, bp_str: str) -> list[Breakpoint]:
        """Parse semicolon-separated breakpoints like 'p11;p11'."""
        bp_parts = bp_str.split(";")
        return [self._parse_breakpoint(bp.strip()) for bp in bp_parts]
