# iscn_authenticator/parser.py
"""Parser for ISCN karyotype strings."""

import re

from iscn_authenticator._abnormality_parsers import _abn, _AbnormalityParsers
from iscn_authenticator._clone_references import resolve_clone_reference
from iscn_authenticator._errors import ParseError
from iscn_authenticator._parser_patterns import (
    DMIN_PATTERN,
    HSR_SIMPLE_PATTERN,
    INCOMPLETE_PATTERN,
    MARKER_PATTERN,
    PSEUDODICENTRIC_PATTERN,
    SEX_CHROMOSOMES_PATTERN,
)
from iscn_authenticator.models import Abnormality, CellLine, KaryotypeAST, Modifiers

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

        original = karyotype
        karyotype = karyotype.strip()
        interphase_prefix = re.match(r"^nuc ish(?:[ \t]+|(?=\())", karyotype)
        if interphase_prefix:
            payload = karyotype[interphase_prefix.end() :]
            return KaryotypeAST(None, "", [], None, Modifiers(ish=payload, interphase=True), raw=original)
        modifiers = Modifiers()
        for marker, attribute in ((".ish ", "ish"), (".arr", "arr"), (".ogm", "ogm")):
            if marker in karyotype:
                karyotype, data = karyotype.split(marker, 1)
                setattr(modifiers, attribute, data)
                break
        mosaic = karyotype.startswith("mos ")
        chimera = karyotype.startswith("chi ")
        if mosaic or chimera:
            karyotype = karyotype[4:]
        ast = self._parse_mosaic(karyotype) if "/" in karyotype else self._parse_single_karyotype(karyotype)
        if mosaic or chimera:
            modifiers.mosaic, modifiers.chimera = mosaic, chimera
        if modifiers != Modifiers():
            if ast.modifiers:
                modifiers.composite = ast.modifiers.composite
            ast.modifiers = modifiers
        ast.raw = original
        return ast

    def _parse_single_karyotype(
        self, karyotype: str, extract_count: bool = False, preceding: list[CellLine] | None = None
    ) -> tuple[KaryotypeAST, int] | KaryotypeAST:
        """Parse a single (non-mosaic) karyotype string."""
        original = karyotype
        count = 0
        composite = False
        if extract_count or "[" in karyotype or "]" in karyotype:
            count_match = re.fullmatch(r"(.+?)(?:\.cp\[([0-9]+)\]|\[(cp)?([0-9]+)\])", karyotype)
            if not count_match and ("[" in karyotype or "]" in karyotype):
                raise ParseError("Cell count must be a positive safe integer in square brackets")
            if count_match:
                karyotype = count_match.group(1)
                composite = count_match[2] is not None or count_match[3] is not None
                digits = (count_match[2] or count_match[4]).lstrip("0") or "0"
                error = (
                    "Composite cell count must be a positive safe integer"
                    if composite
                    else "Cell count must be a positive safe integer in square brackets"
                )
                if len(digits) > 16:
                    raise ParseError(error)
                count = int(digits)
                if not 0 < count <= 9_007_199_254_740_991:
                    raise ParseError(error)

        if "[" in karyotype or "]" in karyotype:
            raise ParseError("Unexpected bracketed count inside the cell line")
        if "," not in karyotype:
            raise ParseError("Missing comma separator between chromosome count and sex chromosomes")

        parts = karyotype.split(",")
        chromosome_count = self._parse_chromosome_count(parts[0])
        token = parts[1].strip()
        sex_chromosomes = token if re.match(r"^(idem|sl|sdl)", token) else self._parse_sex_chromosomes(token)
        abnormalities = self._parse_abnormalities(parts[2:]) if len(parts) > 2 else []

        ast = KaryotypeAST(
            raw=original,
            count=count or None,
            chromosome_count=chromosome_count,
            sex_chromosomes=sex_chromosomes,
            abnormalities=abnormalities,
            cell_lines=None,
            modifiers=Modifiers(composite=count) if composite else None,
        )
        resolve_clone_reference(ast, preceding or [])
        return (ast, count) if extract_count else ast

    def _parse_mosaic(self, karyotype: str) -> KaryotypeAST:
        """Parse a mosaic karyotype with multiple cell lines."""
        groups = karyotype.split("//")
        if len(groups) > 2 or "///" in karyotype:
            raise ParseError("Use one double slash to separate recipient and donor cell lines")
        cell_lines: list[CellLine] = []
        for group_index, group in enumerate(groups):
            if not group.strip() and len(groups) == 2:
                continue
            for line_str in group.split("/"):
                try:
                    ast, count = self._parse_single_karyotype(
                        line_str.strip(),
                        extract_count=True,
                        preceding=[line for line in cell_lines if line.is_donor == (group_index == 1)],
                    )
                except ParseError as error:
                    raise ParseError(f"Cell line {len(cell_lines) + 1}: {error}") from error
                cell_lines.append(
                    CellLine(
                        raw=ast.raw,
                        composite=bool(ast.modifiers and ast.modifiers.composite is not None),
                        reference=ast.reference,
                        changes=ast.changes,
                        chromosome_count=ast.chromosome_count,
                        sex_chromosomes=ast.sex_chromosomes,
                        abnormalities=ast.abnormalities,
                        count=count,
                        is_donor=group_index == 1,
                    )
                )
        if not cell_lines:
            raise ParseError("At least one cell line is required")
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
        # ASCII-only, matching the TypeScript port's /^\d+$/; str.isdigit() also
        # accepts characters such as '²' that int() cannot parse.
        if not (count_str.isascii() and count_str.isdigit()):
            raise ParseError(f"Invalid chromosome count: '{count_str}' is not a number")
        digits = count_str.lstrip("0") or "0"
        if len(digits) > 16:
            raise ParseError("Chromosome count must be a positive safe integer")
        return int(digits)

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
            copy = re.search(r"[x×]([0-9]+)$", stripped)
            token = stripped[: copy.start()] if copy else stripped
            _, clean, uncertain, inheritance = self._strip_modifiers(token)
            abn = self._dispatch_abnormality(clean)
            if copy:
                digits = copy[1].lstrip("0") or "0"
                if len(digits) > 16 or not 0 < int(digits) <= 9_007_199_254_740_991:
                    raise ParseError("Abnormality copy count must be a positive safe integer")
                abn.copy_count = int(digits)
            abn.uncertain = uncertain
            abn.inheritance = inheritance
            abn.raw = stripped
            abnormalities.append(abn)
        return abnormalities

    @staticmethod
    def _strip_modifiers(part: str) -> tuple[str, str, bool, str | None]:
        """Strip leading '?' and trailing inheritance suffix; return (original, clean, uncertain, inheritance)."""
        original = part
        part = part.replace("−", "-")
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
        num_match = re.fullmatch(r"([+-])([0-9]{1,2}|[XY])(c)?", part)
        if num_match:
            return _abn(num_match[1] + (num_match[3] or ""), num_match[2], [], part)

        mar_match = MARKER_PATTERN.match(part)
        if mar_match:
            return self._build_marker(mar_match, part)

        for prefix, method_name in _PREFIX_DISPATCH:
            if part.startswith(prefix):
                return getattr(self, method_name)(part)

        if PSEUDODICENTRIC_PATTERN.match(part):
            return self._parse_pseudodicentric(part)
        outer = re.fullmatch(r"(rec)\(([0-9]{1,2}|[XY])\)(.+)|(upd)\(([0-9]{1,2}|[XY])\)", part)
        if outer:
            return _abn(outer[1] or outer[4], outer[2] or outer[5], [], part)
        complex_outer = re.fullmatch(r"(cpx|cth|cpy)\(([^)]+)\)", part)
        if complex_outer:
            return _abn(complex_outer[1], complex_outer[2], [], part)
        derivative = re.fullmatch(r"([+-]?der)\(([^)]+)\)(.+)", part)
        if derivative:
            return _abn(derivative[1], derivative[2], [], part)
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
