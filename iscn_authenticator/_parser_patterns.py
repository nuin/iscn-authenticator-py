# iscn_authenticator/_parser_patterns.py
"""Compiled regex patterns used by the ISCN karyotype parser.

Kept in a separate module so ``parser.py`` and ``_abnormality_parsers.py``
can share patterns without circular imports.
"""

import re

# Top-level structure
CELL_LINE_COUNT_PATTERN = re.compile(r"^(.+?)\[(\d+)\]$")
SEX_CHROMOSOMES_PATTERN = re.compile(r"^[XYU]+$")

# Breakpoints
BREAKPOINT_PATTERN = re.compile(r"^([pq])(\d+)(?:\.(\d+))?$")
DOUBLE_BREAKPOINT_PATTERN = re.compile(r"^([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)$")
TRIPLE_BREAKPOINT_PATTERN = re.compile(r"^([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)([pq]\d+(?:\.\d+)?)$")

# Numerical abnormalities
NUMERICAL_ABNORMALITY_PATTERN = re.compile(r"^([+-])(\d{1,2}|[XY])$")
MARKER_PATTERN = re.compile(r"^\+(\d*)mar(\d*)$")

# Structural abnormalities (per-type)
DELETION_PATTERN = re.compile(r"^del\((\d{1,2}|[XY])\)\(([^)]+)\)$")
DUPLICATION_PATTERN = re.compile(r"^dup\((\d{1,2}|[XY])\)\(([^)]+)\)$")
INVERSION_PATTERN = re.compile(r"^inv\((\d{1,2}|[XY])\)\(([^)]+)\)$")
TRANSLOCATION_PATTERN = re.compile(r"^t\(([^)]+)\)\(([^)]+)\)$")
ISOCHROMOSOME_SHORT_PATTERN = re.compile(r"^i\((\d{1,2}|[XY])([pq])\)$")
ISOCHROMOSOME_LONG_PATTERN = re.compile(r"^i\((\d{1,2}|[XY])\)\(([^)]+)\)$")
RING_SIMPLE_PATTERN = re.compile(r"^r\((\d{1,2}|[XY])\)$")
RING_BREAKPOINT_PATTERN = re.compile(r"^r\((\d{1,2}|[XY])\)\(([^)]+)\)$")
INSERTION_PATTERN = re.compile(r"^ins\(([^)]+)\)\(([^)]+)\)$")
ADD_PATTERN = re.compile(r"^add\((\d{1,2}|[XY])\)\(([^)]+)\)$")
TRIPLICATION_PATTERN = re.compile(r"^trp\((\d{1,2}|[XY])\)\(([^)]+)\)$")
DICENTRIC_PATTERN = re.compile(r"^dic\(([^)]+)\)\(([^)]+)\)$")
ISODICENTRIC_PATTERN = re.compile(r"^idic\((\d{1,2}|[XY])\)\(([^)]+)\)$")
FRAGILE_SITE_PATTERN = re.compile(r"^fra\((\d{1,2}|[XY])\)\(([^)]+)\)$")
ROBERTSONIAN_PATTERN = re.compile(r"^rob\(([^)]+)\)\(([^)]+)\)$")
QUADRUPLICATION_PATTERN = re.compile(r"^qdp\((\d{1,2}|[XY])\)\(([^)]+)\)$")
DERIVATIVE_PATTERN = re.compile(r"^der\((\d{1,2}|[XY])\)(.+)$")

# Exact-match / "no chromosome" abnormalities
DMIN_PATTERN = re.compile(r"^dmin$")
HSR_SIMPLE_PATTERN = re.compile(r"^hsr$")
HSR_LOCATION_PATTERN = re.compile(r"^hsr\((\d{1,2}|[XY])\)\(([^)]+)\)$")
PSEUDODICENTRIC_PATTERN = re.compile(r"^psu\s*dic\(([^)]+)\)\(([^)]+)\)$")
ACENTRIC_PATTERN = re.compile(r"^ace\((\d{1,2}|[XY])\)\(([^)]+)\)$")
TELOMERIC_ASSOC_PATTERN = re.compile(r"^tas\(([^)]+)\)\(([^)]+)\)$")
FISSION_PATTERN = re.compile(r"^fis\((\d{1,2}|[XY])\)\(([^)]+)\)$")
NEOCENTROMERE_PATTERN = re.compile(r"^neo\((\d{1,2}|[XY])\)\(([^)]+)\)$")
INCOMPLETE_PATTERN = re.compile(r"^inc$")
