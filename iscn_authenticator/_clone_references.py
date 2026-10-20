"""Bounded related-clone resolution (ISCN 2024 section 6.3.4)."""

import re

from iscn_authenticator._errors import ParseError
from iscn_authenticator.models import CellLine, CloneReference, KaryotypeAST


def resolve_clone_reference(ast: KaryotypeAST, preceding: list[CellLine]) -> None:
    token = ast.sex_chromosomes
    if not re.match(r"^(idem|sl|sdl)", token):
        return
    if not re.fullmatch(r"idem|sl|sdl[1-9]\d*|sdl", token):
        raise ParseError(f"Unsupported clone reference '{token}'; spell out this cell line")
    if not preceding:
        raise ParseError(f"Clone reference '{token}' has no preceding reference clone")
    uses_idem = any(line.reference and line.reference.token == "idem" for line in preceding)
    uses_sl = any(line.reference and line.reference.token != "idem" for line in preceding)
    if (token == "idem" and uses_sl) or (token != "idem" and uses_idem):
        raise ParseError("Do not mix idem and sl/sdl clone references")
    sidelines = [line for line in preceding if line.reference and line.reference.token != "idem"]
    reference = preceding[0]
    if token.startswith("sdl"):
        if token == "sdl" and len(sidelines) != 1:
            raise ParseError("Unnumbered sdl clone reference is absent or ambiguous; use sdlN")
        if token != "sdl" and len(token[3:]) > len(str(len(sidelines))):
            raise ParseError(f"Clone reference '{token}' has no preceding sideline")
        index = 0 if token == "sdl" else int(token[3:]) - 1
        if index >= len(sidelines):
            raise ParseError(f"Clone reference '{token}' has no preceding sideline")
        reference = sidelines[index]
    elif any(not line.reference for line in preceding[1:]):
        raise ParseError(f"Clone reference '{token}' is ambiguous after an independent clone")
    ast.reference = CloneReference(token, preceding.index(reference) + 1)
    ast.changes = ast.abnormalities
    ast.sex_chromosomes = reference.sex_chromosomes
    ast.abnormalities = reference.abnormalities + ast.changes
