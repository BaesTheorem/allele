# Roadmap

## The goal is not parity with every genome tool

That target is unreachable, and aiming at it would break the argument this
project rests on.

The constraint is the input, not the effort. A consumer array genotypes a fixed
panel of roughly 600,000 to 1,000,000 positions chosen in advance. It cannot
observe anything that is not on the chip. No amount of code recovers data the
assay never collected, so whole categories of tooling are permanently out of
reach from this input:

| Tool category | Needs | Why an array cannot get there |
|---|---|---|
| Variant calling (GATK, DeepVariant) | sequencing reads | there are no reads, only genotype calls at fixed sites |
| CNV and structural variants | coverage depth, split reads | an array reports a genotype, not a depth |
| Repeat expansions | long reads spanning the repeat | the repeat is not on the chip and could not be spanned if it were |
| Star-allele diplotypes (PharmCAT) | full gene coverage plus phase | arrays cover a few defining SNPs, unphased |
| Novel or private variants | ability to see an unanticipated base | the panel only asks about positions someone chose in advance |

The tool already says this in the one place it comes up. `allele/sources/cpic.py`
opens by stating that it reports a guideline exists for a gene and does not call
star-allele diplotypes, because doing that properly needs phased data and full
allele definitions an array cannot supply, and *"guessing a diplotype would be
worse than not offering one."* That is the whole roadmap in miniature: state the
limit, do not route around it.

Claiming parity with tools that read different data would mean asserting
capabilities the input cannot support. That is the exact failure this project
exists to prevent, and it would be a stranger kind of dishonesty than the one
the plausibility checks guard against, because it would be baked into the
premise instead of into a single finding.

## The goal is parity within the category

The real peers are consumer-array interpreters: Promethease and the array-report
layers of the direct-to-consumer sequencing services. Against those, Allele is
already ahead on two things nothing else in the category does. Every association
carries the study behind it, and every demotion shows its arithmetic instead of
asserting a verdict.

So the roadmap is about closing genuine capability gaps inside that scope, not
about growing into a different tool.

## Ceiling, measured

Worth keeping in view when judging any item below.

- ClinVar GRCh38, restricted to what a chip could ever match (single-base
  substitutions carrying an rsID): **2,664,608 indexed, 155,914 notable**.
- A consumer array overlaps a low-thousands subset of that notable set.
- Indels are structurally excluded. 23andMe writes them as `I` and `D`, which
  carry no sequence, so they can never be matched against a reference allele.

Coverage is not the bottleneck a bigger database would fix. The bottleneck is
the panel.

## Tier 1: worth doing

### Polygenic risk scores

The largest genuine capability gap, and the only item here that is
unambiguously worth doing. Unlike everything in the impossible table above, PRS
works natively on array data: it is a weighted sum over common variants, which
is exactly what a chip measures well.

The PGS Catalog publishes harmonised scoring files for both GRCh37 and GRCh38
as gzipped TSV, with `rsID`, `chr_name`, `chr_position`, `effect_allele`,
`other_allele` and `effect_weight`. The harmonised rsIDs mean this joins the way
every existing source joins.

Two constraints to design around from the start:

- **Licensing.** Scoring files default to the EMBL-EBI terms, but individual
  scores may carry their own, and CC BY-NC-ND appears among them. That is the
  AlphaMissense situation again, so PRS has to be an opt-in source the user
  downloads, never something redistributed in this repo or the web bundle.
- **Ancestry transferability.** PRS transfers across ancestries far worse than
  single-variant findings do, and most published scores were derived in
  European-ancestry cohorts. The GWAS source already labels this problem; a PRS
  feature that did not carry the same warning, prominently, would be worse than
  not shipping one.

A score computed on partial coverage is also not the published score. If the
array covers 60% of a model's variants, that has to be stated as part of the
result, not hidden in a footnote. Closing that coverage gap is the single
legitimate reason to want imputation, and the only consumer of it this roadmap
endorses. See [Imputation, and the firewall it needs](#imputation-and-the-firewall-it-needs).

### Inheritance mode

The principled fix for the class of bug that demoted HFE C282Y. The current
guard is a keyword heuristic on ClinVar's classification text
(`COMMON_BY_DESIGN` in `allele/sources/base.py`), which works for the cases seen
so far and does not generalise.

The real signal is whether the condition is dominant or recessive: a 6% allele
producing 0.36% homozygotes is entirely ordinary for a recessive condition and
impossible for a rare dominant one. That distinction is what the
common-pathogenic rule is actually reaching for.

Verified while writing this: **ClinVar's VCF has no inheritance field.** The
join has to go through `CLNDISDB`, which carries disease cross-references.
Sampling 40,000 records: `CLNDISDB` is present on 95%, of which MedGen appears
on 95% and OMIM on 31%. So MedGen is the join key by coverage, not OMIM, and
this is an integration against an external source rather than one more field to
parse.

### ACMG SF v3.2 secondary findings

A small, authoritative, well-defined answer to "what is actually actionable",
maintained by a body whose whole job is deciding that. For the clinical
category it is a better materiality filter than a magnitude threshold, because
it encodes expert judgement about actionability rather than a proxy for it.

## Tier 2: depth on what is already here

- **Partial star-allele calling** where the array happens to cover a gene's
  defining SNPs, reported with explicit uncertainty. The honest version of this
  says which diplotypes remain possible, not which one you have.
- **Carrier status as a deliberate output.** Being heterozygous for a recessive
  pathogenic variant is a real, useful, non-alarming result, and it currently
  falls out of the pipeline incidentally rather than by design.
- **Molecular consequence.** `MC=` is present on 100% of sampled ClinVar
  records and is not parsed at all today. Cheap to add, and it tells the reader
  whether a variant is missense, nonsense or synonymous.
- **`CLNSIGCONF`**, which holds the breakdown behind a conflicting
  classification. Conflicts are currently filtered out by `is_notable`. Showing
  the split would fit the project's stance better than dropping it.

## Imputation, and the firewall it needs

Shipping as an installable app removes the packaging objection to imputation.
It does not remove the statistical one, and the statistical one is the reason to
be careful.

First, a correction worth stating plainly: the browser build was never what
blocked this. The CLI already installs locally and runs offline, so it could
have carried a heavy dependency all along. What an installable app actually buys
is permission to drop "standard library only", which is a different and smaller
thing than it sounds.

### Imputation is least accurate exactly where clinical findings live

Imputation accuracy is a function of allele frequency, and the relationship is
steep. From the Haplotype Reference Consortium paper: imputing from an Illumina
1M array at **MAF 0.1%**, HRC reaches **r² = 0.64** and 1000 Genomes Phase 3
reaches **r² = 0.36**. The standard quality threshold used in that same analysis
is **r² ≥ 0.5**.

So with the only panel that is freely redistributable, a variant at 0.1%
frequency already fails the accepted quality bar. Now compare that against this
codebase's own arithmetic: `COHORT_DETECTION_LIMIT` in
`allele/sources/base.py` is `1e-5`, and the whole cohort-absence rule is about
variants rarer than roughly 1 in 100,000. That is two orders of magnitude below
the point where free-panel imputation has already stopped working.

Running ClinVar annotation over imputed genotypes would therefore manufacture
precisely the false positives the entire plausibility architecture exists to
suppress, and it would manufacture them with a confident provenance trail
attached. The better panels (HRC, TOPMed) are gated behind access applications,
so they cannot ship with an offline tool anyway.

### The rule

**Imputed calls never feed clinical annotation.** ClinVar, CPIC and
AlphaMissense read genotyped calls only.

Imputation exists to serve polygenic scores, which is the one consumer of it
that the accuracy curve supports: PRS is a weighted sum dominated by common
variants, and common variants are where imputation is reliable. That makes the
dependency worth carrying for exactly one feature.

Any imputed call that reaches a report at all carries its `r²`, and the
threshold is shown rather than assumed. A genotype the model guessed is a
different kind of claim from a genotype the chip measured, and the report should
never present the two as the same.

### Suggested shape: three tiers, not one replacement

Collapsing the browser build into an app would trade a real asset for
convenience. Keep them layered instead:

| Tier | Install | Dependencies | Scope |
|---|---|---|---|
| Web | none, drop a file on a page | none | shipped subset, findings only |
| CLI | `pip install -e .` | standard library, `pysam` for gnomAD | full databases, current behaviour |
| App | installer, downloads panels | numpy plus Beagle or Minimac | adds imputation and PRS |

The web build is the zero-install artifact that makes this shareable, it already
works, and `tests/test_parity.py` now keeps it honest. Losing it to gain a
feature that only the top tier needs would be a bad trade.

One thing worth weighing before adding the third tier: "standard library only"
is not just an aesthetic. Every dependency is code with read access to the most
sensitive file a person owns, in a tool whose pitch is that the genome never
leaves the machine. That argues for the heavy tier being opt-in and clearly
bounded rather than the default install, in the same way gnomAD and
AlphaMissense already are.

## Non-goals

- Anything from the impossible table. Not deferred, not out of reach for now:
  structurally unavailable from this input.
- Diagnostic framing of any kind. See Scope in the README.
- Redistributing a database whose licence does not permit it. Currently
  AlphaMissense and SNPedia, and prospectively parts of the PGS Catalog.
- Resolving disagreements between sources. Showing the disagreement is the
  feature.
- Imputed genotypes reaching ClinVar, CPIC or AlphaMissense. Imputation is
  unreliable in the frequency range those sources care about, so this stays a
  hard boundary rather than a threshold to tune.
- Retiring the zero-install web build to fund the app tier. They are layers, not
  alternatives.
