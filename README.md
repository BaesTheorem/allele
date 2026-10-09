# Allele

Read a consumer DNA export and see what the public databases actually say about
it. Runs entirely offline. Your genome never leaves the machine.

**Web version: [baestheorem.github.io/allele](https://baestheorem.github.io/allele)** —
runs in the browser tab, no upload, no server.

```bash
pip install -e .
allele db update --build 37 && allele db update --build 38   # both assemblies
allele report ~/Downloads/genome.txt                          # writes HTML
allele report genome.txt -f json -o r.json                    # or JSON / VCF
allele report genome.txt --panel my_snps.txt                  # restrict to a panel
allele diff old.json new.json                                 # what ClinVar changed
allele db status | allele db path | allele db clean
```

## Inputs

Detected from file content, not the extension, because vendors ship
tab-delimited data in files named `.csv`. Plain, gzipped and zipped all work.

| Input | Notes |
|---|---|
| 23andMe | 4 columns, concatenated genotype |
| AncestryDNA | 5 columns, separate alleles; chromosomes 23/24/25/26 map to X/Y/X/MT |
| MyHeritage | Quoted CSV |
| Family Tree DNA | Quoted CSV, same shape |
| Living DNA | Tab-delimited despite the `.csv` name |
| VCF | Single sample, or the first of a multi-sample file with a warning |
| Promethease report | Decodes the embedded payload plus SNPedia's curation |

## Sources

| Source | Records | License |
|---|---|---|
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | 2,664,744 rsID-bearing SNVs | Public domain |
| [GWAS Catalog](https://www.ebi.ac.uk/gwas/) | 847,508 associations, 341,074 variants | Free, attribution |
| [CPIC](https://cpicpgx.org) | 979 variants, 19 genes, levels A/B | Free |
| [gnomAD](https://gnomad.broadinstitute.org) v4.1 exomes | frequency for notable positions | ODbL 1.0 |
| [AlphaMissense](https://github.com/google-deepmind/alphamissense) | 226,998 predictions | CC BY-NC-SA, not redistributed |
| [SNPedia](https://www.snpedia.com) | 43,779 variants | CC BY-NC-SA, read from your own report |

Two of these carry non-commercial, share-alike terms, so neither is ever
redistributed here. AlphaMissense is downloaded by the user on request. SNPedia
is read only from a Promethease report already on your disk, which embeds the
wiki's grading of every genotype at the variants it covers.

gnomAD and AlphaMissense are opt-in and GRCh38-only:

```bash
pip install -e ".[gnomad]"              # gnomAD needs pysam; nothing else does
allele db update --only gnomad          # targeted range requests, ~1 hour
allele db update --only alphamissense   # 0.6 GB download, indexes in a minute
```

Both are keyed by GRCh38 coordinate and neither will answer for a GRCh37 file.
A position means nothing without the assembly it was measured on, and a GRCh37
coordinate lands on an unrelated GRCh38 variant often enough to matter, so both
sources stay silent rather than attach a frequency or a prediction belonging to
something you do not have.

gnomAD is an enricher rather than a source of findings. It supplies the
population frequency the plausibility checks need; "this variant is common" is
not a result on its own. AlphaMissense only speaks where ClinVar is silent,
since a prediction adds nothing beside an expert-panel classification, and is
always labelled a prediction rather than a finding.

Where sources disagree, the disagreement is shown rather than resolved. ClinVar
reclassifying something a year-old SNPedia grading calls harmless is the useful
part, not a conflict to paper over.

## Every association cites its study

A p-value with no provenance is unreadable: it gives you no way to tell a
754-person study from 2007 apart from a 1.6-million-person meta-analysis from
2025, and no way to go read either. So each trait association carries the paper
behind it.

- **PubMed link** with first author, journal and year — 100% coverage
- **Discovery cohort size** — 100% coverage
- **Ancestry of that cohort** — 85% coverage

Two flags come straight out of that. 60% of catalog associations were
discovered in European-ancestry cohorts only, which is labelled, because effect
sizes frequently do not transfer. Discovery cohorts under 5,000 are labelled
too, because small studies overstate effects.

ClinVar findings carry their review status as a star rating, which is the
single best signal for how much agreement sits behind a classification.

## Why most of this code is about being wrong

A consumer array reporting a rare pathogenic variant is usually wrong. That is
the central problem, and on real data it dominated everything else: before these
checks, the top of the ranked output was almost entirely genotyping artifacts,
including variants SNPedia annotates in plain words as *"Ancestry miscall in all
likelihood"* and *"99% of Ancestry users are told they have this genotype"*.

Three checks demote them, each showing its arithmetic rather than asserting a
verdict.

**Hardy-Weinberg.** An allele at frequency *f* appears homozygous in about *f²*
of people. When that expected rate falls far below the array's own error rate,
an error is the better explanation.

**Absent from every cohort.** A pathogenic variant with no frequency in gnomAD,
1000 Genomes, ExAC or ESP is missing because it is rarer than those cohorts can
resolve. That bounds *f* below roughly 1 in 100,000, so a homozygous or
hemizygous call is almost certainly an artifact. This is the fallback; where
gnomAD has a measured frequency it is used instead, and the explanation says
which source the number came from.

**Common pathogenic.** A variant classified pathogenic yet carried by more than
5% of people cannot cause a rare dominant condition. That is usually a
reference-allele or annotation error.

Demoted findings are still shown, under their own heading with the reasoning.
Hiding them would be its own kind of dishonesty.

### Where these rules were wrong about you

Each check is an argument about a *rare, highly penetrant, dominant* condition,
and each one had to learn where that premise fails.

Frequency is the load-bearing number, and there is no such thing as one
frequency. HFE C282Y sits near 1% globally and near 6% in Northern Europeans,
where it is the diagnostic genotype for hereditary haemochromatosis in roughly 1
in 200 people. Judged against the global average, a true homozygote failed
Hardy-Weinberg and was filed as a genotyping error. So the frequency used is the
highest across gnomAD's ancestry groups, not the average: the question is
whether *this reader* could carry the genotype, and their ancestry is unknown.

Being common is also not always evidence against a classification. ClinVar marks
some variants low penetrance, a risk factor or a drug-response marker, and those
are *expected* to be common. C282Y is classified
`Pathogenic/Pathogenic, low penetrance|risk factor`, so the common-pathogenic
rule now stands down where ClinVar's own wording predicts a common allele, and
the explanation names recessive inheritance as the other reading rather than
asserting an artifact.

A measured frequency of exactly zero used to mean "no opinion", because the
cohort-absence fallback tested for a missing value rather than a falsy one.
Indexing gnomAD therefore switched the check off for thousands of the variants
it was written for, which made the tool safer when its enrichment source was
broken than when it worked.

## The other things that go quietly wrong

**Strand.** SNPedia reports a large minority of variants on the minus strand:
19,428 of 56,087 in one real report, about one in three. ClinVar and the GWAS
Catalog are plus-strand, so the canonical model is plus-strand and conversions
are explicit. In Promethease reports the field that drives the flip is
`flipped`, not `orientation` — the two disagree on 2,314 variants in that same
report. A/T and C/G heterozygotes read identically on both strands and are
flagged rather than guessed.

**Ploidy.** 23andMe writes a Y or mitochondrial call as `G`; AncestryDNA writes
`G G`. Both are haploid, so the single allele is honest and the doubled form is
a layout artifact. Without normalizing, one person's two exports disagree at
every Y and MT site. X is left alone, because collapsing it means inferring sex.

**Build.** Headers lie. Both ClinVar assemblies are indexed, so any variant
whose GRCh37 and GRCh38 positions differ becomes a probe: a sample of the file's
own coordinates is checked against both, and the coordinates win over the header
when they disagree. On the first real file tested, a Promethease report
declaring `reference: 37` turned out to carry GRCh38 coordinates, 86 probes to
0. Matching is by rsID so annotations were unaffected, but any position-based
logic would have been off by megabases.

**Carriage.** A variant is only a finding if you carry the allele. A pathogenic
variant you do not have is not a result.

**Materiality.** SNPedia grades nearly every variant it knows, most of them 0
for "boring". A floor requiring a ClinVar entry with assertion criteria,
magnitude 2+, or a carried genome-wide significant association cuts 48,576 raw
annotations to about 4,800.

## Web version

The browser build ships a subset so it can run with no server: ClinVar entries
with at least one review star (both assemblies, with positions so it can verify
build the same way the CLI does), the strongest few associations per variant,
and all of CPIC. About 10 MB gzipped.

It omits gnomAD and AlphaMissense. AlphaMissense cannot be redistributed at
all, and gnomAD frequency for the full notable set would dominate the download.
The command line version therefore reports more, though not strictly more per
finding: without gnomAD the browser falls back to ClinVar's sparser frequency
fields, so the two can reach different verdicts on the same variant.

Both implementations follow the same pipeline, and `tests/test_parity.py`
enforces it rather than asserting it. Identical inputs go through both engines
and the verdicts are compared: the plausibility rules across every combination
of zygosity, frequency and classification, plus scoring, zygosity, genotype and
chromosome normalization, ploidy collapsing and strand ambiguity. The
comparison is rule-level, because the CLI annotates against the full database
and the browser against a shipped subset, so different *findings* are expected
and different *rules* are not.

That test was written after the claim turned out to be false. It found
pharmacogenomic findings scoring 45 in the CLI and 0 in the browser, and a
combined `Pathogenic/Likely pathogenic` classification that the JavaScript
ranked correctly and the Python sent to the unrecognized default. It skips when
Node is unavailable.

## Development

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[dev]"
.venv/bin/python -m pytest
.venv/bin/python scripts/build_web_bundle.py   # regenerate docs/db
```

Standard library only, with one exception: gnomAD reads a remote VCF over
htslib range requests, so `allele db update --only gnomad` needs
`pip install -e ".[gnomad]"`. Nothing else imports it, and every other command
works without it.

The parity tests additionally need Node, and skip without it. Vendor fixtures
are synthetic; real exports, generated reports and downloaded databases are all
gitignored.

## Scope

Not a diagnostic. Consumer arrays test a fixed set of positions, miss most of
the genome, and produce false positives at rare sites. Nothing here should
drive a medical decision. Confirm anything that matters with a clinical-grade
test and a genetic counselor.

### The limits are structural, not a backlog

An array genotypes a panel of positions somebody chose in advance, roughly
600,000 to 1,000,000 of them. It cannot see anything else, so a whole class of
questions is not answered slowly here, it is unanswerable from this input. No
version of this tool will call a novel variant, a copy-number change or a repeat
expansion, because there are no sequencing reads to call them from. It will not
report a pharmacogenomic star-allele diplotype either, which needs full gene
coverage and phase, and is why the CPIC findings say a guideline exists for the
gene rather than claiming your genotype changes a dose.

Even the annotation ceiling is lower than the database sizes suggest. ClinVar
has 155,914 notable variants among the single-base substitutions a chip could
match at all; a consumer array overlaps a few thousand of them. Indels are out
by construction: 23andMe writes them as `I` and `D`, which carry no sequence, so
they can never be checked against a reference allele.

Reading "0 clinical findings" as "nothing to worry about" is therefore the
wrong inference. It means nothing was found among the small, pre-selected set of
positions this file happens to contain.

[ROADMAP.md](ROADMAP.md) covers what is worth building given that ceiling, and
what is permanently outside it.

## License

MIT for the code. Databases keep their own terms, listed above.
