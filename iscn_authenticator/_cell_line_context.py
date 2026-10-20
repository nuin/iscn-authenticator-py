"""Build a group context without leaking a composite count to other groups."""

from dataclasses import replace

from iscn_authenticator.models import CellLine, KaryotypeAST, Modifiers


def cell_line_context(ast: KaryotypeAST, line: CellLine) -> KaryotypeAST:
    modifiers = replace(ast.modifiers or Modifiers(), composite=line.count if line.composite else None)
    return replace(
        ast,
        chromosome_count=line.chromosome_count,
        sex_chromosomes=line.sex_chromosomes,
        abnormalities=line.abnormalities,
        cell_lines=None,
        raw=line.raw,
        count=line.count,
        reference=line.reference,
        changes=line.changes,
        modifiers=modifiers if modifiers != Modifiers() else None,
    )
