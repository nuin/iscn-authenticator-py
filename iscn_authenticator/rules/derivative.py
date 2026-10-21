"""Checks of flat named components, without derivative path reconstruction."""

import re
from dataclasses import replace

from iscn_authenticator._errors import ParseError
from iscn_authenticator.parser import KaryotypeParser
from iscn_authenticator.rules.routine import structural_reference_rule


def _balanced_parentheses(text: str) -> bool:
    depth = 0
    for char in text:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


def validate_derivative_components(ast, abn, rules) -> list[str]:
    if abn.type not in {"der", "+der", "-der"}:
        return []
    if not _balanced_parentheses(abn.raw):
        return [f"Unbalanced parentheses in {abn.raw}"]
    tail = abn.raw[abn.raw.index(")") + 1 :]
    components = re.finditer(
        r"\b(del|dup|inv|t|ins|add|trp|qdp|dic|idic|fra|rob|hsr|ace|tas|fis|neo|i|r)"
        r"\(([^()]*)\)(?:\(([^()]*)\))?",
        tail,
    )
    errors = _direct_breakpoints(ast, abn, tail)
    for match in components:
        try:
            if match[3] is None:
                component = replace(abn, type=match[1], chromosome=match[2], raw=match[0])
                errors.extend(structural_reference_rule.validate(ast, component))
            else:
                component = KaryotypeParser()._parse_abnormalities([match[0]])[0]
                errors.extend(error for rule in rules for error in rule.validate(ast, component))
        except ParseError as error:
            errors.append(str(error))
    return [f"Derivative component in {abn.raw}: {error}" for error in errors]


def _direct_breakpoints(ast, abn, tail) -> list[str]:
    match = re.match(r"\(([pq0-9.;?]+)\)", tail)
    if not match:
        return []
    try:
        breakpoints = KaryotypeParser()._parse_multiple_breakpoints(match[1])
        if len(breakpoints) != len(abn.chromosome.split(";")):
            return [f"Derivative chromosome and breakpoints counts differ in {abn.raw}"]
        return structural_reference_rule.validate(ast, replace(abn, breakpoints=breakpoints))
    except ParseError as error:
        return [str(error)]
