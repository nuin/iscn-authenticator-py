"""Bounded related-clone resolution (ISCN 2024 section 6.3.4)."""

import re

from iscn_authenticator._clone_arithmetic import expand_clone
from iscn_authenticator._errors import ParseError
from iscn_authenticator.models import CellLine, CloneReference, KaryotypeAST


def resolve_clone_reference(ast: KaryotypeAST, preceding: list[CellLine]) -> None:
    token = ast.sex_chromosomes
    if not re.match(r"^(idem|sl|sdl)", token):
        return
    match = re.fullmatch(r"(idem|sl[1-9][0-9]*|sl|sdl[1-9][0-9]*|sdl)(?:[x×]([0-9]+))?", token)
    if not match:
        raise ParseError(f"Unsupported clone reference '{token}'; spell out this cell line")
    base = match[1]
    digits = (match[2] or "1").lstrip("0") or "0"
    if len(digits) > 16 or not 0 < int(digits) <= 9_007_199_254_740_991:
        raise ParseError("Clone multiplier must be a positive safe integer")
    multiplier = int(digits)
    if not preceding:
        raise ParseError(f"Clone reference '{token}' has no preceding reference clone")
    uses_idem = any(line.reference and line.reference.token.startswith("idem") for line in preceding)
    uses_sl = any(line.reference and not line.reference.token.startswith("idem") for line in preceding)
    if (base == "idem" and uses_sl) or (base != "idem" and uses_idem):
        raise ParseError("Do not mix idem and sl/sdl clone references")
    sidelines = [line for line in preceding if line.reference and not line.reference.token.startswith("idem")]
    stemlines = [line for line in preceding if not line.reference]
    reference = preceding[0]
    if base.startswith("sdl"):
        if len(stemlines) > 1:
            raise ParseError(
                "Sideline references with multiple stemlines require explicit restatement; spell out this cell line"
            )
        if base == "sdl" and len(sidelines) != 1:
            raise ParseError("Unnumbered sdl clone reference is absent or ambiguous; use sdlN")
        if base != "sdl" and len(base[3:]) > len(str(len(sidelines))):
            raise ParseError(f"Clone reference '{token}' has no preceding sideline")
        index = 0 if base == "sdl" else int(base[3:]) - 1
        if index >= len(sidelines):
            raise ParseError(f"Clone reference '{token}' has no preceding sideline")
        reference = sidelines[index]
    elif base.startswith("sl") and base != "sl":
        if len(base[2:]) > len(str(len(stemlines))) or int(base[2:]) > len(stemlines):
            raise ParseError(f"Clone reference '{token}' has no preceding stemline")
        reference = stemlines[int(base[2:]) - 1]
    elif base == "sl" and len(stemlines) > 1:
        raise ParseError(f"Clone reference '{token}' is ambiguous after an independent clone")
    if reference.composite or (ast.modifiers and ast.modifiers.composite is not None):
        raise ParseError("Composite groups cannot use clone references; spell out each group")
    ast.reference = CloneReference(token, preceding.index(reference) + 1)
    ast.changes = ast.abnormalities
    adjusted = expand_clone(ast, reference, multiplier)
    if multiplier != 1:
        ast.reference.multiplier = multiplier
    if adjusted:
        ast.reference.adjusted = True
