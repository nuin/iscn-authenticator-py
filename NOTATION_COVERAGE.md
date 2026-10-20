# Notation coverage

The validator checks a bounded subset of ISCN 2024. A result with `valid: true` means the implemented checks found no errors. It does not establish full ISCN compliance, biological clonality, clinical significance, or completeness of a laboratory report. Use the coverage result alongside validity.

## Coverage API

TypeScript exports `getValidationCoverage(parsed)`; Python exports `get_validation_coverage(parsed)`. Both return `level`, `notes`, `reference`, `rule_version`, and `band_reference`. Levels are `supported`, `partial`, and `unparsed`. `supported` describes coverage of the implemented routine checks and is independent of whether those checks passed. A parsed invalid karyotype can therefore have supported coverage.

The reference is `ISCN 2024`; the rule version is `2026-10-19`; the band reference is `UCSC hg38 cytoBand`. Applications must present partial or unparsed results for review. The existing `valid`, `errors`, and `parsed` validation contract is unchanged.

## Checked and deferred notation

| Notation | Implemented checks | Remaining limits |
|---|---|---|
| Routine single line, e.g. `46,XX[20]` | Numeric chromosome count, sex field, supported abnormalities, optional positive metaphase count | A bracketed count does not make the input a mosaic |
| Multiple cell lines | Each line receives its own rules; each requires a positive metaphase count; normal lines follow abnormal lines | Ordering multiple abnormal clones requires specimen context and further tie-break rules |
| Recipient/donor groups | One `//` boundary; groups and line order preserved | Donor identity and transplant context require review |
| Numerical autosomal changes | Chromosomes 1–22; conservative total check where the diploid assumption is reliable | No inferred total for ranges, composite descriptions, unknown or uncertain changes, multiplicities, or count-changing structural notation |
| Acquired sex chromosome changes | X/Y identifiers, syntax, explicit loss-of-X context including a remaining Y | Total coherence is deferred: acquired changes and the reported sex complement cannot be naively added together |
| Polyploid/hypohaploid counts | Positive safe integers and ascending ranges; counts are not restricted to 23–92 | Ploidy and full numerical balance require review; positive arithmetic alone does not establish biological plausibility |
| Common structural abnormalities | Existing breakpoint-number/arm checks, chromosome identifiers, mapped band checks | This does not check all event-specific biological or nomenclature constraints |
| Insertions | One chromosome with three concatenated breakpoints, or two chromosomes with a site and a two-breakpoint segment separated by `;` | Other insertion layouts are explicitly rejected as unsupported; no unmatched layout bypasses band checks |
| Related clones | `count,idem,...`, `count,sl,...`, and `count,sdlN,...`; unnumbered `sdl` only with one preceding sideline | Multiple stemlines, multipliers on references, and ambiguous/absent references are rejected for explicit restatement |
| Composite counts | Final `[cpN]`, with legacy `.cp[N]` accepted; positive safe count | Per-abnormality frequencies and mixed composite/noncomposite line grammar are not fully supported |
| `der`, `rec`, `cpx`, `cth`, `cpy`, `upd` | Outer form and chromosome identifiers | Nested rearrangements, derivation, phase, and complex structure are partial coverage |
| Unknown abnormality tokens | Original text retained | These produce partial coverage even if no implemented rule reports an error |
| FISH | Standalone `nuc ish` and appended `.ish` payload retained | Full FISH grammar, signal interpretation, probe validity and clinical cutoffs are not validated by core |
| Array and OGM | Appended `.arr` / `.ogm` payload retained | Full assay grammar, coordinates, assemblies and interpretation are not validated; standalone assay descriptions may be unparsed |

The laboratory review FISH builder checks its separate, bounded structured copy-number input. That does not broaden core coverage of arbitrary FISH expressions. CNV evidence assessment is also separate from notation validity.

## Numerical totals

The automatic total uses 44 autosomes plus the written sex complement and explicit autosomal gains/losses. It applies only to numeric counts from 40 through 55, one to four X/Y characters, and a restricted set of count-preserving autosomal structural types (`del`, `dup`, `inv`, `t`, `add`, `trp`, `qdp`, `i`, `fra`, `ins`). Any numerical X/Y change, structural sex-chromosome change, unknown event, multiplicity, uncertainty, composite, or other structural type defers the total. This is a deliberately bounded check, not a general ploidy reconstruction algorithm.

For example, `46,XX,+21` fails the implemented total check; `47,XX,+21[20]` passes it. A polyploid input such as `184,XXXXXXXX[4]` retains its count and receives partial coverage. A reported total outside the diploid-check interval receives no inferred balance verdict. Zero, negative, fractional, or unsafe chromosome and metaphase counts are rejected. Cross-runtime integers are limited to 9,007,199,254,740,991.

## Cytoband reference

The bundled catalogue is the UCSC hg38 `cytoBand` dataset, reused by the ideogram. It is not a complete ISCN atlas and does not prove the existence of every finer ISCN subband. A catalogue match includes an exact band or a parent band supported by catalogue children. A designation finer than a known catalogue band is retained with a coverage note. A region-only breakpoint such as `p1?` retains `band: null` and uncertainty; it is never converted to centromeric `p10`.

Literal `p10` and `q10` are recognized centromeric designations. A band that matches neither a catalogue band nor its hierarchy is reported as absent from this particular catalogue; users should verify reference disagreement against the authoritative ISCN atlas. Malformed band digits such as `p101` are rejected rather than normalized into a different band. The validator checks the correct chromosome for each breakpoint in supported insertions and multi-chromosome events.

Catalogue source: [UCSC hg38 cytoBand download](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz). Python bundles an identical copy of the TypeScript catalogue, without a runtime download or dependency.

## Related-clone provenance

`idem` refers to the first clone in its recipient/donor group, not the immediately previous line. `sl` refers to that group's first stemline. Sidelines are numbered in their written order: `sdl1` references the first preceding resolved sideline, `sdl2` the second. Unnumbered `sdl` requires exactly one preceding sideline. Mixing `idem` with `sl`/`sdl`, referring across the donor boundary, or guessing a stemline after an independent clone is rejected. The implementation does not infer biological clone identity from similarity.

The AST keeps the original complete input in `raw`. Single-line metaphase counts use `count` and leave `cell_lines: null`. Each explicit cell line keeps `raw`; a resolved related line also carries `reference: { token, cell_line }`, where `cell_line` is one-based within its recipient/donor group, and `changes`, containing the written additional abnormalities. Its `abnormalities` list contains the resolved inherited plus additional findings. Each abnormality retains its literal `raw` text, including accepted Unicode minus and multiplication signs. This expansion represents the written reference and does not establish that the clones are biologically related.

## References and verification

The implementation was checked against local ISCN 2024 sections 4.5.3, 5.3.1.2, 5.4, 6.3.4–6.3.7 and these publisher corrections:

- [March 2025 erratum, DOI 10.1159/000543493](https://karger.com/cgr/article-pdf/165/2/99/4350726/000543493.pdf): includes corrected composite-example bands.
- [July 2025 erratum, DOI 10.1159/000546684](https://karger.com/cgr/article/165/6/369/929293/Erratum): includes a corrected polyploid range and assay examples.
- [December 2025 online / 2026 issue erratum, DOI 10.1159/000549238](https://karger.com/cgr/article/166/1/64/939149/Erratum-ISCN-2024-An-International-System-for): includes corrected chromosome totals and derivative descriptions.

The errata review does not imply complete implementation of every corrected example or assay. Synthetic regression cases cover the audited failures, neighbouring valid inputs, ambiguous clone references, acquired sex losses, polyploid counts, band precision, insertion mapping, uncertainty and oversized numbers. The shared immutable fixture corpus is consumed by both runtimes; Python also runs Hypothesis properties. See the release guide for the exact committed dependency and corpus pins.
