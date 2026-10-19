"""Checks for explicit cell lines (ISCN 2024 section 4.5.3)."""

from iscn_authenticator.models import CellLine, KaryotypeAST
from iscn_authenticator.rules.base import Rule


def is_normal_cell_line(line: CellLine) -> bool:
    return line.chromosome_count == 46 and line.sex_chromosomes in {"XX", "XY"} and not line.abnormalities


def validate_cell_line_counts(ast: KaryotypeAST, target) -> list[str]:
    return [
        f"Cell line {index}: Cell count must be a positive safe integer in square brackets"
        for index, line in enumerate(ast.cell_lines or [], 1)
        if type(line.count) is not int or not 0 < line.count <= 9_007_199_254_740_991
    ]


def validate_cell_line_order(ast: KaryotypeAST, target) -> list[str]:
    errors: list[str] = []
    for donor in (False, True):
        previous_normal = None
        for index, line in enumerate(ast.cell_lines or [], 1):
            if line.is_donor != donor:
                continue
            if not is_normal_cell_line(line):
                if previous_normal:
                    errors.append(
                        f"Cell line {index}: The normal cell line must be listed last within each recipient "
                        "or donor group (ISCN 2024 section 4.5.3)."
                    )
                continue
            if previous_normal and (
                previous_normal.count < line.count
                or (
                    previous_normal.count == line.count
                    and previous_normal.sex_chromosomes == "XY"
                    and line.sex_chromosomes == "XX"
                )
            ):
                errors.append(
                    f"Cell line {index}: Normal cell lines must be ordered by decreasing cell count, "
                    "with XX before XY when counts are equal (ISCN 2024 section 4.5.3)."
                )
            previous_normal = line
    return errors


ALL_CELL_LINE_RULES = [
    Rule(
        id="CELL_LINE_COUNTS",
        category="cell_lines",
        description="Every cell line must have a positive integer metaphase count",
        validate=validate_cell_line_counts,
    ),
    Rule(
        id="CELL_LINE_ORDER",
        category="cell_lines",
        description="Normal lines follow abnormal lines within recipient and donor groups",
        validate=validate_cell_line_order,
    ),
]
