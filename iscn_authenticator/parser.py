# iscn_authenticator/parser.py
"""Parser for ISCN karyotype strings."""

from iscn_authenticator._abnormality_parsers import _abn, _AbnormalityParsers
from iscn_authenticator._errors import ParseError
from iscn_authenticator._parser_patterns import (
    CELL_LINE_COUNT_PATTERN,
    DERIVATIVE_PATTERN,
    DMIN_PATTERN,
    HSR_SIMPLE_PATTERN,
    INCOMPLETE_PATTERN,
    MARKER_PATTERN,
    NUMERICAL_ABNORMALITY_PATTERN,
    PSEUDODICENTRIC_PATTERN,
    SEX_CHROMOSOMES_PATTERN,
)
from iscn_authenticator.models import Abnormality, CellLine, KaryotypeAST

__all__ = ["KaryotypeParser", "ParseError"]


_INHERITANCE_SUFFIXES = (("mat", 3), ("pat", 3), ("dn", 2))

# Ordered prefix → parser-method dispatch. More-specific prefixes come before
# their shorter relatives so they win the lookup (e.g. "idic(" before any "i").
_PREFIX_DISPATCH: tuple[tuple[str, str], ...] = (
    ("hsr(", "_parse_hsr_location"),
    ("idic(", "_parse_isodicentric"),
    ("ins(", "_parse_insertion"),
    ("inv(", "_parse_inversion"),
    ("del(", "_parse_deletion"),
    ("dup(", "_parse_duplication"),
    ("dic(", "_parse_dicentric"),
    ("add(", "_parse_add"),
    ("fra(", "_parse_fragile_site"),
    ("trp(", "_parse_triplication"),
    ("qdp(", "_parse_quadruplication"),
    ("rob(", "_parse_robertsonian"),
    ("ace(", "_parse_acentric"),
    ("tas(", "_parse_telomeric_assoc"),
    ("fis(", "_parse_fission"),
    ("neo(", "_parse_neocentromere"),
    ("t(", "_parse_translocation"),
    ("r(", "_parse_ring"),
    ("i(", "_parse_isochromosome"),
)


class KaryotypeParser(_AbnormalityParsers):
    """Parses ISCN karyotype strings into an AST."""

    def parse(self, karyotype: str) -> KaryotypeAST:
        """Parse a karyotype string into an AST."""
        if not karyotype or not karyotype.strip():
            raise ParseError("Karyotype string is empty")

        karyotype = karyotype.strip()
        if "/" in karyotype:
            return self._parse_mosaic(karyotype)
        return self._parse_single_karyotype(karyotype)

    def _parse_single_karyotype(
        self, karyotype: str, extract_count: bool = False
    ) -> tuple[KaryotypeAST, int] | KaryotypeAST:
        """Parse a single (non-mosaic) karyotype string."""
        count = 0
        if extract_count:
            count_match = CELL_LINE_COUNT_PATTERN.match(karyotype)
            if count_match:
                karyotype = count_match.group(1)
                count = int(count_match.group(2))

        if "," not in karyotype:
            raise ParseError("Missing comma separator between chromosome count and sex chromosomes")

        parts = karyotype.split(",")
        chromosome_count = self._parse_chromosome_count(parts[0])
        sex_chromosomes = self._parse_sex_chromosomes(parts[1])
        abnormalities = self._parse_abnormalities(parts[2:]) if len(parts) > 2 else []

        ast = KaryotypeAST(
            chromosome_count=chromosome_count,
            sex_chromosomes=sex_chromosomes,
            abnormalities=abnormalities,
            cell_lines=None,
            modifiers=None,
        )
        return (ast, count) if extract_count else ast

    def _parse_mosaic(self, karyotype: str) -> KaryotypeAST:
        """Parse a mosaic karyotype with multiple cell lines."""
        cell_lines = []
        for line_str in karyotype.split("/"):
            ast, count = self._parse_single_karyotype(line_str.strip(), extract_count=True)
            cell_lines.append(
                CellLine(
                    chromosome_count=ast.chromosome_count,
                    sex_chromosomes=ast.sex_chromosomes,
                    abnormalities=ast.abnormalities,
                    count=count,
                    is_donor=False,
                )
            )

        first = cell_lines[0]
        return KaryotypeAST(
            chromosome_count=first.chromosome_count,
            sex_chromosomes=first.sex_chromosomes,
            abnormalities=first.abnormalities,
            cell_lines=cell_lines,
            modifiers=None,
        )

    def _parse_chromosome_count(self, count_str: str) -> int | str:
        """Parse chromosome count (number or '~'-separated range)."""
        count_str = count_str.strip()
        if "~" in count_str:
            return count_str
        if not count_str.isdigit():
            raise ParseError(f"Invalid chromosome count: '{count_str}' is not a number")
        return int(count_str)

    def _parse_sex_chromosomes(self, sex_str: str) -> str:
        """Parse sex chromosome designation."""
        sex_str = sex_str.strip()
        if not SEX_CHROMOSOMES_PATTERN.match(sex_str):
            raise ParseError(f"Invalid sex chromosomes: '{sex_str}' must contain only X, Y, or U")
        return sex_str

    def _parse_abnormalities(self, parts: list[str]) -> list[Abnormality]:
        """Parse a list of abnormality strings into Abnormality objects."""
        abnormalities: list[Abnormality] = []
        for raw_part in parts:
            stripped = raw_part.strip()
            if not stripped:
                continue
            original, clean, uncertain, inheritance = self._strip_modifiers(stripped)
            abn = self._dispatch_abnormality(clean)
            abn.uncertain = uncertain
            abn.inheritance = inheritance
            abn.raw = original
            abnormalities.append(abn)
        return abnormalities

    @staticmethod
    def _strip_modifiers(part: str) -> tuple[str, str, bool, str | None]:
        """Strip leading '?' and trailing inheritance suffix; return (original, clean, uncertain, inheritance)."""
        original = part
        uncertain = False
        if part.startswith("?"):
            uncertain = True
            part = part[1:]
        inheritance: str | None = None
        for suffix, length in _INHERITANCE_SUFFIXES:
            if part.endswith(suffix):
                inheritance = suffix
                part = part[:-length]
                break
        return original, part, uncertain, inheritance

    def _dispatch_abnormality(self, part: str) -> Abnormality:
        """Map a (modifier-stripped) abnormality string to an Abnormality node."""
        num_match = NUMERICAL_ABNORMALITY_PATTERN.match(part)
        if num_match:
            return _abn(num_match.group(1), num_match.group(2), [], part)

        mar_match = MARKER_PATTERN.match(part)
        if mar_match:
            return self._build_marker(mar_match, part)

        for prefix, method_name in _PREFIX_DISPATCH:
            if part.startswith(prefix):
                return getattr(self, method_name)(part)

        if PSEUDODICENTRIC_PATTERN.match(part):
            return self._parse_pseudodicentric(part)
        if DERIVATIVE_PATTERN.match(part):
            return self._parse_derivative(part)
        if DMIN_PATTERN.match(part):
            return _abn("dmin", "", [], part)
        if HSR_SIMPLE_PATTERN.match(part):
            return _abn("hsr", "", [], part)
        if INCOMPLETE_PATTERN.match(part):
            return _abn("inc", "", [], part)

        return _abn("unknown", "", [], part)

    @staticmethod
    def _build_marker(match, part: str) -> Abnormality:
        count_prefix, marker_suffix = match.group(1), match.group(2)
        chromosome = "mar" + marker_suffix if marker_suffix else "mar"
        abn = _abn("+mar", chromosome, [], part)
        abn.copy_count = int(count_prefix) if count_prefix else None
        return abn
