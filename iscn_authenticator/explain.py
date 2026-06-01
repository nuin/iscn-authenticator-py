"""Human-readable explanations of karyotype AST nodes.

Provides template-based explanations only. Curated explanation lookup
(against a JSON database of clinical-condition mappings) is deferred to a
future release.
"""

from .models import Abnormality, ExplainResult, KaryotypeAST, KaryotypeNode

ABNORMALITY_NAMES = {
    "+": "Gain",
    "-": "Loss",
    "del": "Deletion",
    "dup": "Duplication",
    "inv": "Inversion",
    "t": "Translocation",
    "i": "Isochromosome",
    "r": "Ring chromosome",
    "ins": "Insertion",
    "add": "Additional material",
    "trp": "Triplication",
    "dic": "Dicentric chromosome",
    "idic": "Isodicentric chromosome",
    "fra": "Fragile site",
    "rob": "Robertsonian translocation",
    "mar": "Marker chromosome",
}


def format_breakpoints(node: Abnormality) -> str:
    """Formats breakpoints into a human-readable string."""
    if not node.breakpoints:
        return ""
    bps = []
    for bp in node.breakpoints:
        s = bp.arm + str(bp.region or "") + str(bp.band or "")
        if bp.subband:
            s += "." + bp.subband
        bps.append(s)
    return f" at {', '.join(bps)}"


def generate_template_explanation(node: KaryotypeNode) -> ExplainResult:
    """Generates a deterministic, mechanical description of an AST node."""
    summary = ""
    detail = ""

    if isinstance(node, Abnormality):
        type_name = ABNORMALITY_NAMES.get(node.type, node.type)
        bp_text = format_breakpoints(node)

        if node.type in ("+", "-"):
            summary = f"{type_name} of chromosome {node.chromosome}."
            detail = f"The karyotype indicates a {type_name.lower()} of an entire chromosome {node.chromosome}."
        elif node.type == "mar":
            summary = "Marker chromosome."
            detail = "An unidentified extra structurally abnormal chromosome (ESAC) is present."
        else:
            summary = f"{type_name} on chromosome {node.chromosome}{bp_text}."
            detail = f"A {type_name.lower()} was identified on chromosome {node.chromosome}{bp_text}."

        if node.inheritance:
            inh_map = {
                "mat": "maternally inherited",
                "pat": "paternally inherited",
                "dn": "de novo (not inherited)",
            }
            inh_text = inh_map.get(node.inheritance, f"inherited ({node.inheritance})")
            summary += f" ({node.inheritance})"
            detail += f" This abnormality is {inh_text}."

    elif isinstance(node, KaryotypeAST):
        count = node.chromosome_count
        sex = node.sex_chromosomes
        abncount = len(node.abnormalities)

        summary = f"{count},{sex} karyotype with {abncount} abnormalities."
        detail = f"This is a {sex} karyotype with a total chromosome count of {count}. "
        if abncount == 0:
            detail += "No structural or numerical abnormalities were detected."
        else:
            detail += f"There are {abncount} abnormality/abnormalities described."

    return ExplainResult(summary=summary, detail=detail, citation=None, refs={}, confidence="template")


def explain(node: KaryotypeNode) -> ExplainResult:
    """Explains a karyotype AST node in human-readable terms.

    Currently returns a template-based explanation. A curated-explanation
    lookup (with confidence="curated") is planned for a future release.
    """
    return generate_template_explanation(node)
