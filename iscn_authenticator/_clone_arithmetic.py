"""Bounded clone multiplication and unique autosomal structural subtraction."""

import re
from dataclasses import replace

from iscn_authenticator._errors import ParseError
from iscn_authenticator.models import Abnormality, CellLine, KaryotypeAST

_MAX_SAFE = 9_007_199_254_740_991
_NEUTRAL = {"del", "dup", "inv", "t", "add", "trp", "qdp", "i", "fra", "ins"}


def _copies(abn: Abnormality, count: int) -> Abnormality:
    if not 0 < count <= _MAX_SAFE:
        raise ParseError("Expanded abnormality copy count exceeds the safe integer limit")
    return replace(abn, copy_count=None if count == 1 else count)


def _subtract_structure(inherited: list[Abnormality], change: Abnormality) -> Abnormality:
    from iscn_authenticator.parser import KaryotypeParser

    token = change.raw.replace("−", "-")
    match = re.fullmatch(
        r"-(del|dup|inv|add|trp|qdp|i|r|idic)\(([1-9]|1[0-9]|2[0-2])\)"
        r"(?:\(([^()]*)\))?(?:[x×][0-9]+)?",
        token,
    )
    if not match or change.uncertain:
        raise ParseError(f"Unsupported structural clone subtraction '{change.raw}'; spell out this cell line")
    specified = (
        KaryotypeParser()._parse_abnormalities([f"{match[1]}({match[2]})({match[3]})"])[0]
        if match[3] is not None
        else None
    )
    candidates = [
        abn
        for abn in inherited
        if abn.type == match[1]
        and abn.chromosome == match[2]
        and (specified is None or abn.breakpoints == specified.breakpoints)
    ]
    if len(candidates) != 1:
        reason = "ambiguous" if candidates else "absent from the reference clone"
        raise ParseError(f"Structural clone subtraction '{change.raw}' is {reason}; spell out this cell line")
    selected = candidates[0]
    if selected.uncertain or any(bp.uncertain for bp in selected.breakpoints):
        raise ParseError(f"Cannot subtract uncertain inherited structure '{change.raw}'; spell out this cell line")
    remaining = (selected.copy_count or 1) - (change.copy_count or 1)
    if remaining < 0:
        raise ParseError(f"Structural clone subtraction '{change.raw}' exceeds inherited copies")
    index = inherited.index(selected)
    if remaining:
        inherited[index] = _copies(selected, remaining)
    else:
        inherited.pop(index)
    return replace(change, type="-", chromosome=match[2], breakpoints=[])


def _append_numerical(findings: list[Abnormality], change: Abnormality) -> bool:
    remaining = change.copy_count or 1
    adjusted = False
    if change.type in {"+", "-"} and re.fullmatch(r"[1-9]|1[0-9]|2[0-2]", change.chromosome) and not change.uncertain:
        for candidate in findings[:]:
            if (
                candidate.type != ("-" if change.type == "+" else "+")
                or candidate.chromosome != change.chromosome
                or candidate.uncertain
                or candidate.inheritance != change.inheritance
            ):
                continue
            removed = min(remaining, candidate.copy_count or 1)
            count = (candidate.copy_count or 1) - removed
            index = findings.index(candidate)
            if count:
                findings[index] = _copies(candidate, count)
            else:
                findings.pop(index)
            remaining -= removed
            adjusted = True
            if not remaining:
                break
    if remaining:
        findings.append(_copies(change, remaining))
    return adjusted


def _count_change(abn: Abnormality) -> int | None:
    if abn.uncertain or any(bp.uncertain for bp in abn.breakpoints):
        return None
    if abn.type in {"+", "-"} and re.fullmatch(r"[1-9]|1[0-9]|2[0-2]", abn.chromosome):
        return (1 if abn.type == "+" else -1) * (abn.copy_count or 1)
    return 0 if abn.type in _NEUTRAL and not re.search(r"[XY]", abn.chromosome) else None


def expand_clone(ast: KaryotypeAST, reference: CellLine, multiplier: int) -> bool:
    sex = reference.sex_chromosomes
    if len(sex) * multiplier > 4096:
        raise ParseError("Clone sex complement exceeds the 4096-character expansion limit; spell out this cell line")
    ast.sex_chromosomes = "U" if sex == "U" else "".join(sorted(sex * multiplier))
    inherited = [_copies(abn, (abn.copy_count or 1) * multiplier) for abn in reference.abnormalities]
    additions = []
    adjusted = bool(reference.reference and (reference.reference.adjusted or reference.reference.multiplier))
    delta = 0
    for original in ast.changes or []:
        structural_loss = bool(re.match(r"\??[-−][a-z]", original.raw))
        change = _subtract_structure(inherited, original) if structural_loss else original
        adjusted = adjusted or structural_loss
        additions.append(change)
        value = _count_change(change)
        delta = None if delta is None or value is None else delta + value
    for change in additions:
        cancelled = _append_numerical(inherited, change)
        adjusted = adjusted or cancelled
    ast.abnormalities = inherited
    if isinstance(reference.chromosome_count, int) and (multiplier != 1 or adjusted):
        base = reference.chromosome_count * multiplier
        expected = None if delta is None else base + delta
        if not 0 < base <= _MAX_SAFE or (expected is not None and not 0 < expected <= _MAX_SAFE):
            raise ParseError("Expanded clone chromosome count exceeds the positive safe integer limit")
        if expected is not None and isinstance(ast.chromosome_count, int) and expected != ast.chromosome_count:
            raise ParseError(
                f"Chromosome count {ast.chromosome_count} does not match the expected total {expected} "
                "from the reference clone and written changes"
            )
    return adjusted
