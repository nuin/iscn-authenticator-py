"""Conservative count and bundled UCSC hg38 cytoband checks."""

import json
import re
from pathlib import Path

from iscn_authenticator.models import Abnormality, KaryotypeAST
from iscn_authenticator.rules.base import Rule

BANDS = json.loads(Path(__file__).resolve().parents[1].joinpath("bands.json").read_text())
STABLE_COUNT_TYPES = {"del", "dup", "inv", "t", "add", "trp", "qdp", "i", "fra", "ins"}


def expected_diploid_count(ast: KaryotypeAST) -> int | None:
    if ast.reference and (ast.reference.multiplier or ast.reference.adjusted):
        return None
    count = ast.chromosome_count
    if (
        not isinstance(count, int)
        or not 40 <= count <= 55
        or not re.fullmatch(r"[XY]{1,4}", ast.sex_chromosomes)
        or (ast.modifiers and ast.modifiers.composite is not None)
    ):
        return None
    expected = 44 + len(ast.sex_chromosomes)
    for abn in ast.abnormalities:
        if abn.uncertain or abn.copy_count is not None:
            return None
        if abn.type in {"+", "-", "+c", "-c"}:
            if abn.chromosome in {"X", "Y"}:
                return None
            expected += 1 if abn.type.startswith("+") else -1
        elif abn.type not in STABLE_COUNT_TYPES or re.search(r"[XY]", abn.chromosome):
            return None
    return expected


def breakpoint_name(bp) -> str:
    return f"{bp.arm}{bp.region}{bp.band if bp.band is not None else ''}" + (f".{bp.subband}" if bp.subband else "")


def parent_band(parent: str, child: str) -> bool:
    return child == parent or child.startswith(parent + ".") or ("." in parent and child.startswith(parent))


def band_status(chromosome, bp) -> str:
    name = breakpoint_name(bp)
    if name in {"p10", "q10"}:
        return "known"
    bands = BANDS.get(chromosome, [])
    if bp.band is None:
        return "finer" if any(band["name"].startswith(name) for band in bands) else "absent"
    if any(parent_band(name, band["name"]) for band in bands):
        return "known"
    if any(parent_band(band["name"], name) for band in bands):
        return "finer"
    return "absent"


def chromosome_breakpoints(abn: Abnormality):
    chromosomes = abn.chromosome.split(";")
    if len(chromosomes) == 1:
        return [(chromosomes[0], bp) for bp in abn.breakpoints]
    if abn.type == "ins" and len(chromosomes) == 2 and len(abn.breakpoints) == 3:
        return [(chromosomes[0 if index == 0 else 1], bp) for index, bp in enumerate(abn.breakpoints)]
    if len(chromosomes) == len(abn.breakpoints):
        return list(zip(chromosomes, abn.breakpoints, strict=True))
    return []


def validate_structural_reference(_, abn: Abnormality) -> list[str]:
    if not abn.chromosome or abn.chromosome.startswith("mar") or abn.type in {"unknown", "+", "-", "+c", "-c"}:
        return []
    invalid = [
        chromosome for chromosome in abn.chromosome.split(";") if not re.fullmatch(r"[1-9]|1\d|2[0-2]|X|Y", chromosome)
    ]
    if invalid:
        return [f"Invalid chromosome '{chromosome}' in {abn.raw}. Must be 1-22, X, or Y" for chromosome in invalid]
    return [
        f"Band '{chromosome}{breakpoint_name(bp)}' in {abn.raw} is absent from the bundled UCSC hg38 cytoband catalogue"
        for chromosome, bp in chromosome_breakpoints(abn)
        if band_status(chromosome, bp) == "absent"
    ]


def validate_clone_position(_, abn) -> list[str]:
    if re.fullmatch(r"idem|sl\d*|sdl\d*", abn.raw):
        return ["Related-clone shorthand must replace the sex chromosome field and reference a preceding clone"]
    return []


def validate_numerical_total(ast, _) -> list[str]:
    expected = expected_diploid_count(ast)
    if expected is not None and expected != ast.chromosome_count:
        return [
            f"Chromosome count {ast.chromosome_count} does not match the expected total {expected} "
            "from the sex complement and numerical changes"
        ]
    return []


structural_reference_rule = Rule(
    "ABN_REFERENCE", "abnormality", "Bundled chromosome and band reference", validate_structural_reference
)
misplaced_clone_rule = Rule(
    "ABN_CLONE_POSITION", "abnormality", "Related-clone shorthand position", validate_clone_position
)
numerical_total_rule = Rule(
    "CHR_NUMERICAL_TOTAL", "coherence", "Reliable diploid numerical total", validate_numerical_total
)
