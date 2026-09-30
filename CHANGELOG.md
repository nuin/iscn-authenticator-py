# Changelog

All notable changes to `iscn-authenticator` are recorded here. The format
roughly follows [Keep a Changelog](https://keepachangelog.com/); versions
are [SemVer](https://semver.org/) and on `0.x` a breaking change bumps
the minor.

Sections per release: **Added**, **Changed**, **Fixed**, **Removed**.

---

## Unreleased

### Added
- Hypothesis property tests (`tests/test_properties.py`) that generate
  valid and invalid karyotypes; `hypothesis` is now in the `test` extra.

### Fixed
- `validate_karyotype` no longer raises `ValueError` on chromosome counts
  made of non-ASCII digit characters such as `²`; counts are now ASCII-only,
  matching the TypeScript port.

## 0.2.1 — 2026-06

### Changed
- Internal refactor of `parser.py` and `rules/abnormality.py`: extracted
  regex patterns, per-type abnormality parsers, and validator
  implementations into private modules; replaced the monolithic
  abnormality dispatch with a small prefix-driven table. Public API
  (`KaryotypeParser`, `ParseError`, rule instances) is unchanged.
- Project homepage URL updated to
  `https://bioinformat.org/projects/iscn-authenticator`.

---

## 0.2.0 — 2026-04

First release tracked in this changelog. First release published to PyPI.

### Added
- ISCN 2024 grammar coverage: numerical aberrations (`+`/`-`), deletions
  (`del`), duplications (`dup`), translocations (`t`), inversions
  (`inv`), insertions (`ins`), isochromosomes (`i`/`idic`),
  derivatives (`der`/`dic`), rings (`r`), Robertsonian translocations
  (`rob`), triplications (`trp`), marker chromosomes (`mar`),
  uncertainty (`?`), inheritance suffixes (`mat`/`pat`/`dn`),
  mosaicism (cell lines split on `/`).
- AST + rule-engine architecture (`KaryotypeParser` →
  `KaryotypeAST` → `RuleEngine`); rules split into chromosome-level
  and abnormality-level lists.
- `validate_karyotype(s) -> ValidationResult` (`{ valid, errors, parsed }`)
  and `is_valid_karyotype(s) -> bool` convenience wrapper.
- Shared fixture corpus at `fixtures/validity.json` exercising both this
  library and the TypeScript port.

### Deferred to a later release
- Curated explanation lookup (clinical-condition database mapping karyotype
  signatures to authoritative summaries). The template-based explanation
  via `generate_template_explanation` and `explain` is included; the
  curated lookup will return in a future minor release.

### Notes
- Zero runtime dependencies; standard library only.
- Supports Python 3.10–3.13.
