# iscn-authenticator

ISCN 2024 karyotype validation — parser, rule engine, and AST types. Pure Python, zero runtime dependencies, supports Python 3.10–3.13.

> Validates karyotype strings against the 2024 International System for Human Cytogenomic Nomenclature (ISCN) and returns a structured AST alongside any rule violations.

## Install

```bash
pip install iscn-authenticator
```

## Quick start

```python
from iscn_authenticator import is_valid_karyotype, validate_karyotype

is_valid_karyotype("46,XX")                              # True
is_valid_karyotype("47,XY,+21")                          # True
is_valid_karyotype("46,XX,t(9;22)(q34;q11.2)")           # True
is_valid_karyotype("47,XY,+21[8]/46,XY[12]")             # True (mosaic)
is_valid_karyotype("garbage")                            # False

result = validate_karyotype("46,XX,del(5)(q13q33)")
# result.valid   -> True
# result.errors  -> []
# result.parsed  -> KaryotypeAST(...)
```

## API surface

- `validate_karyotype(s: str) -> ValidationResult` — full validation, returns `{valid, errors, parsed}`.
- `is_valid_karyotype(s: str) -> bool` — boolean convenience wrapper.
- `KaryotypeAST`, `Abnormality`, `Breakpoint`, `CellLine` — dataclasses describing the parsed structure.
- `ParseError` — raised by the parser on malformed syntax; `validate_karyotype` catches it and surfaces it through `errors`.

## Supported nomenclature

ISCN 2024 grammar covered:

- Numerical aberrations: `+21`, `-X`, ...
- Structural rearrangements: `del`, `dup`, `inv`, `t` (translocation), `ins` (insertion).
- Specialized rearrangements: `i` (isochromosome), `idic` (isodicentric), `der` (derivative), `dic` (dicentric), `r` (ring), `rob` (Robertsonian), `trp` (triplication).
- Marker chromosomes: `mar`.
- Uncertainty (`?`), inheritance suffixes (`mat`, `pat`, `dn`), mosaicism (cell lines split on `/`).

## ISCN 2024 reference

Implementation tracks the 2024 edition of the International System for Human Cytogenomic Nomenclature. See the [ISCN publication](https://karger.com/books/book/277/ISCN-2024).

## Cross-implementation parity

A TypeScript port (`@iscn/core`) lives in the [iscn-authenticator monorepo](https://github.com/nuin/iscn-authenticator) under `packages/core/`, and consumes the same fixture corpus shipped in this repo at `fixtures/validity.json`. Parity is enforced via CI in both repos. The Hypothesis property tests in `tests/test_properties.py` generate their own cases and never write to the shared corpus.

## Contributing

```bash
git clone https://github.com/nuin/iscn-authenticator-py.git
cd iscn-authenticator-py
pip install pre-commit ruff
pip install -e ".[test]"   # hypothesis, for tests/test_properties.py
pre-commit install
python3 -m unittest discover tests
```

## License

MIT — see `LICENSE`.
