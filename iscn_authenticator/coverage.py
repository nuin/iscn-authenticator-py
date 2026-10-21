"""Coverage metadata for the bounded ISCN checks; independent of validity."""

from iscn_authenticator._cell_line_context import cell_line_context
from iscn_authenticator.models import KaryotypeAST
from iscn_authenticator.rules.cell_lines import is_normal_cell_line
from iscn_authenticator.rules.routine import band_status, chromosome_breakpoints, expected_diploid_count


def get_validation_coverage(ast: KaryotypeAST | None) -> dict:
    notes = {}
    if ast is None:
        notes["unparsed"] = "The notation could not be parsed; no complete structure was checked."
    else:
        modifiers = ast.modifiers
        if modifiers and any(getattr(modifiers, field) is not None for field in ("ish", "arr", "ogm")):
            notes["assay"] = (
                "FISH, array and OGM payloads are preserved but their complete assay grammar is not validated."
            )
        if (modifiers and modifiers.composite is not None) or any(line.composite for line in ast.cell_lines or []):
            notes["composite"] = (
                "Composite counts are checked; per-abnormality frequencies and chromosome totals require review."
            )
        lines = [cell_line_context(ast, line) for line in ast.cell_lines] if ast.cell_lines else [ast]
        for line in lines:
            if line.reference and (line.reference.multiplier or line.reference.adjusted):
                notes["clone_resolution"] = (
                    "Clone multiplication or subtraction was expanded within bounded rules; "
                    "review the original notation and resolved findings together."
                )
            if line.chromosome_count is not None and expected_diploid_count(line) is None:
                notes["count_scope"] = (
                    "A reliable diploid numerical total cannot be inferred for this notation; "
                    "review ploidy, structural changes and chromosome balance."
                )
            for abn in line.abnormalities:
                if abn.type == "unknown":
                    notes["unknown"] = "Unrecognized abnormality text is preserved and requires manual review."
                if abn.type in {"der", "+der", "-der", "rec", "cpx", "cth", "cpy", "upd"}:
                    notes["complex"] = (
                        "Recognizable derivative components receive chromosome, band and breakpoint checks; "
                        "complete derivative paths, recombinant and complex structures require review."
                    )
                if abn.uncertain or any(bp.uncertain for bp in abn.breakpoints):
                    notes["uncertain"] = "Uncertain abnormalities or breakpoints require review."
                if any(band_status(chromosome, bp) == "finer" for chromosome, bp in chromosome_breakpoints(abn)):
                    notes["band_resolution"] = (
                        "A region-only or finer subband designation requires review "
                        "beyond the bundled catalogue resolution."
                    )
        if any(
            sum(not is_normal_cell_line(line) for line in ast.cell_lines or [] if line.is_donor == donor) > 1
            for donor in (False, True)
        ):
            notes["clone_order"] = (
                "Ordering of multiple abnormal clones requires specimen context and further ISCN tie-break rules."
            )
        if any(line.is_donor for line in ast.cell_lines or []):
            notes["donor"] = (
                "Recipient and donor groups are preserved; donor identity and transplant context require review."
            )
    return {
        "level": "unparsed" if ast is None else "partial" if notes else "supported",
        "notes": [{"code": code, "message": message} for code, message in notes.items()],
        "reference": "ISCN 2024",
        "rule_version": "2026-10-20.2",
        "band_reference": "UCSC hg38 cytoBand",
    }
