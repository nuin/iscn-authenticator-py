# iscn_authenticator/models.py
from dataclasses import dataclass


@dataclass
class ExplainResult:
    """Result of the Explain module."""

    summary: str  # One sentence summary
    detail: str  # Detailed description
    citation: dict | None  # ISCN spec citation: {"section": "...", "page": ...}
    refs: dict[str, list[str]]  # External references: {"omim": [...]}
    confidence: str  # "template" | "curated" | "none"


@dataclass
class Breakpoint:
    """Represents a chromosomal breakpoint like p11.2 or q34."""

    arm: str  # "p", "q", "cen", "ter"
    region: int | None
    band: int | None
    subband: str | None  # "1", "11", etc.
    uncertain: bool  # Has ? in designation


@dataclass
class Abnormality:
    """Represents a chromosomal abnormality."""

    type: str  # "+", "-", "del", "dup", "inv", "t", etc.
    chromosome: str  # "5" or "9;22" for translocations
    breakpoints: list[Breakpoint]
    inheritance: str | None  # "mat", "pat", "dn", etc.
    uncertain: bool  # Has ? marker
    copy_count: int | None  # For x2, x3 notation
    raw: str  # Original string


@dataclass
class Modifiers:
    """Karyotype-level modifiers."""

    mosaic: bool = False  # mos
    chimera: bool = False  # chi
    constitutional: bool = False  # c suffix
    incomplete: bool = False  # inc


@dataclass
class CellLine:
    """Represents a cell line in mosaic/chimera notation."""

    chromosome_count: int
    sex_chromosomes: str
    abnormalities: list[Abnormality]
    count: int  # Number in brackets [10]
    is_donor: bool = False  # For chimera: after //


@dataclass
class KaryotypeAST:
    """Abstract syntax tree for a parsed karyotype."""

    chromosome_count: int | str  # int or range "45~48"
    sex_chromosomes: str  # "XX", "XY", "X", etc.
    abnormalities: list[Abnormality]
    cell_lines: list[CellLine] | None
    modifiers: Modifiers | None


@dataclass
class ValidationResult:
    """Result of karyotype validation."""

    valid: bool
    errors: list[str]
    parsed: KaryotypeAST | None
    explanation: ExplainResult | None = None


KaryotypeNode = KaryotypeAST | Abnormality
