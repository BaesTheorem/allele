# genome-report

Read a consumer DNA export and produce an annotated health and traits report.
Entirely offline: your genome never leaves the machine.

> **Status: early.** Input parsing is built and tested. The annotation engine
> and the report renderer are next. See [Roadmap](#roadmap).

Working name. Not published yet.

## Supported inputs

Detected from file content, not from the extension, because vendors are
careless about it. Plain, gzipped and zipped files all work.

| Input | Notes |
|---|---|
| 23andMe | 4 columns, concatenated genotype. Build from header. |
| AncestryDNA | 5 columns, separate allele columns. Chromosomes 23/24/25/26 map to X/Y/X/MT. |
| MyHeritage | Quoted CSV. |
| Family Tree DNA | Quoted CSV, same shape as MyHeritage. |
| Living DNA | Tab-delimited despite the `.csv` name. |
| VCF | Single sample, or the first sample of a multi-sample file with a warning. |
| Promethease report | Decodes the embedded genotype payload plus SNPedia's curation. |

## Why the Promethease path exists

Promethease is closed source, but the report it generates is self-contained
HTML with the data separated from the code. Genotypes sit in base64'd,
zlib-compressed JSON. Decoding that recovers the calls, SNPedia's curation of
them, and, uniquely among these inputs, the strand each SNP was reported on.

A report also embeds how SNPedia grades *every* genotype at the SNPs it covers,
not only the one you carry: in a real report, 43,779 SNPs and 60,821 genotype
gradings. So a report you already own can annotate a raw export from the same
person.

That table is SNPedia content under CC BY-NC-SA. Reading it out of a file on
your own disk is fine. Redistributing it is not, so it is extracted at runtime
and never cached into this package.

## The things that quietly go wrong

These are the failure modes the code is built around, all of them silent if
unhandled.

**Strand.** SNPedia reports a large minority of SNPs on the minus strand. In a
real report, 19,428 of 56,087 oriented SNPs, about one in three. Join on rsID
without checking orientation and a third of genotypes read backwards. A/T and
C/G heterozygotes are worse: both strands look valid, so they cannot be
oriented by allele identity at all and are flagged rather than guessed.

**Ploidy.** 23andMe writes a Y or mitochondrial call as `G`. AncestryDNA's
two-column layout writes the same call as `G G`. Y and MT are haploid, so the
single allele is honest and the doubled form is a layout artifact. Without
normalizing, one person's two files disagree at every Y and MT site. X is
deliberately left alone, because collapsing it means inferring sex and being
wrong about that is worse.

**Build.** Positions are meaningless without knowing the assembly. When a
header declares no build, this refuses to assume one and says so, rather than
defaulting to 37 and being quietly wrong on a build 38 export.

**Indels.** 23andMe encodes insertions and deletions as `I` and `D`, which
carry no sequence and can never be matched to a reference allele. They are
marked unusable instead of being compared to something they cannot equal.

**No-calls.** `--`, `00`, `NN` and `.` all mean the chip produced nothing.
Dropped, not treated as data.

## Usage

```python
from genome_report import parse

sample = parse("~/Downloads/genome.txt")
print(sample.source_format, sample.build, len(sample))
print(sample.stats())
for warning in sample.warnings:
    print("warning:", warning)
```

## Development

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[dev]"
.venv/bin/python -m pytest
```

Vendor fixtures under `fixtures/` are synthetic. Real exports are personal data
and are gitignored, along with generated reports and downloaded databases.

Standard library only. No pandas: these files stream, and loading one into a
dataframe on a laptop is a waste of memory.

## Roadmap

1. ~~Input parsing, normalization, vendor detection~~ done
2. Annotation from local databases: ClinVar for clinical significance, GWAS
   Catalog for trait associations. Downloaded once, indexed into SQLite.
3. Report renderer: self-contained HTML, no network at view time.
4. Confidence handling. A consumer chip calling a rare pathogenic variant is
   more often a genotyping error than a real finding, which is why 23andMe
   limits its own health reports and why Promethease shows a false-positive
   warning on rare high-magnitude hits. Anything pathogenic has to be presented
   with that caveat attached, never as a result.

## Scope and limits

Not a diagnostic. Consumer genotyping chips are not clinical sequencing: they
test a fixed set of positions, miss most of the genome, and produce false
positives on rare variants. Nothing here should be used to make a medical
decision. Confirm anything that matters with a clinical-grade test and a
genetic counselor.

## License

MIT for the code. Annotation databases keep their own licenses: ClinVar is
public domain, the GWAS Catalog is free to use with attribution, and SNPedia
content is CC BY-NC-SA and is never redistributed here.
