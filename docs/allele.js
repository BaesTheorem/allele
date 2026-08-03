/* Allele -- browser engine.
 *
 * This is a port of the Python pipeline, kept deliberately parallel to it:
 * parse -> annotate -> plausibility -> materiality -> render. It runs entirely
 * in the page. There is no upload endpoint, no fetch of user data, no
 * analytics. The only network requests are for the static annotation bundles
 * that ship alongside this file.
 */
'use strict';

/* ------------------------------------------------------------ genotypes -- */

const CHROMS = new Set(['1','2','3','4','5','6','7','8','9','10','11','12','13',
  '14','15','16','17','18','19','20','21','22','X','Y','MT']);
const CHROM_ALIASES = {'23':'X','24':'Y','25':'X','26':'MT','M':'MT','MITO':'MT','XY':'X'};
const NO_CALL = new Set(['--','-','00','0','NN','N','','..','.']);
const NUCS = new Set(['A','C','G','T']);
const INDELS = new Set(['I','D']);
const COMPLEMENTS = {A:'T',T:'A',C:'G',G:'C'};
const HAPLOID = new Set(['Y','MT']);

function normalizeChromosome(raw) {
  if (raw == null) return null;
  let v = String(raw).trim().toUpperCase();
  if (v.startsWith('CHR')) v = v.slice(3);
  v = CHROM_ALIASES[v] || v;
  return CHROMS.has(v) ? v : null;
}

function normalizeGenotype(raw) {
  if (raw == null) return '';
  let v = String(raw).trim().toUpperCase().replace(/[();|/,\s]/g, '');
  if (NO_CALL.has(v)) return '';
  for (const ch of v) if (!NUCS.has(ch) && !INDELS.has(ch)) return '';
  return v.split('').sort().join('');
}

// Y and MT are haploid; AncestryDNA's two-column layout doubles them.
function normalizePloidy(chrom, genotype) {
  if (HAPLOID.has(chrom) && genotype.length === 2 && genotype[0] === genotype[1]) {
    return genotype[0];
  }
  return genotype;
}

function complement(g) {
  return g.split('').map(c => COMPLEMENTS[c] || c).sort().join('');
}

// A/T and C/G heterozygotes read the same on either strand.
function isStrandAmbiguous(g) {
  if (g.length !== 2 || g[0] === g[1]) return false;
  const s = new Set(g);
  return (s.has('A') && s.has('T')) || (s.has('C') && s.has('G'));
}

function isIndel(g) { for (const c of g) if (INDELS.has(c)) return true; return false; }
function usable(call) { return !!call.g && !isIndel(call.g); }

/* -------------------------------------------------------------- parsing -- */

const BUILD_RE = /(?:build|grch)\s*([0-9]{2})/i;
const BUILD_MAP = {'36':36,'37':37,'38':38,'19':37,'18':36};

function sniffBuild(text) {
  const m = text.match(BUILD_RE);
  return m ? (BUILD_MAP[m[1]] || null) : null;
}

function detectFormat(header, firstLine) {
  const h = header.toLowerCase();
  if (h.includes('23andme')) return '23andme';
  if (h.includes('ancestrydna') || h.includes('ancestry.com')) return 'ancestrydna';
  if (h.includes('myheritage')) return 'myheritage';
  if (h.includes('living dna') || h.includes('livingdna')) return 'livingdna';
  if (h.includes('ftdna') || h.includes('family tree dna')) return 'ftdna';
  const f = (firstLine || '').toLowerCase();
  if (f.startsWith('rsid,chromosome,position,result')) return 'ftdna';
  if (f.startsWith('rsid\tchromosome\tposition\tallele1')) return 'ancestrydna';
  return 'generic';
}

function splitLine(line) {
  if (line.includes('\t')) return line.split('\t').map(s => s.trim().replace(/^"|"$/g, ''));
  const out = []; let cur = ''; let quoted = false;
  for (const ch of line) {
    if (ch === '"') { quoted = !quoted; continue; }
    if (ch === ',' && !quoted) { out.push(cur.trim()); cur = ''; continue; }
    cur += ch;
  }
  out.push(cur.trim());
  return out;
}

function parseChipText(text) {
  const lines = text.split(/\r?\n/);
  const comments = [];
  let i = 0;
  for (; i < lines.length && i < 80; i++) {
    if (lines[i].startsWith('#')) comments.push(lines[i]); else break;
  }
  const header = comments.join(' ');
  const format = detectFormat(header, lines[i]);
  const build = sniffBuild(header);
  const warnings = [];
  if (build === null) {
    warnings.push('No genome build declared in the header. Positions cannot be ' +
      'verified, so annotation falls back to rsID matching only.');
  }

  const calls = [];
  let malformed = 0;
  for (; i < lines.length; i++) {
    const line = lines[i];
    if (!line || line.startsWith('#')) continue;
    const f = splitLine(line);
    if (f.length < 4) { if (line.trim()) malformed++; continue; }
    const first = f[0].toLowerCase();
    if (first === 'rsid' || first === 'rs id' || first === 'snp') continue;
    const chrom = normalizeChromosome(f[1]);
    const raw = f.length === 4 ? f[3] : f[3] + f[4];
    const pos = parseInt(f[2], 10);
    calls.push({
      rsid: f[0],
      chrom,
      pos: Number.isNaN(pos) ? null : pos,
      g: normalizePloidy(chrom, normalizeGenotype(raw)),
      build,
    });
  }
  if (malformed) warnings.push(`${malformed.toLocaleString()} malformed line(s) skipped`);
  return { calls, format, build, warnings, orientation: {}, flipped: {}, curation: {} };
}

function parseVcf(text) {
  const calls = [];
  const warnings = [];
  let build = null;
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith('##')) {
      if (build === null && (line.includes('reference=') || line.includes('assembly='))) {
        build = sniffBuild(line);
      }
      continue;
    }
    if (line.startsWith('#CHROM')) {
      const cols = line.split('\t');
      if (cols.length > 10) {
        warnings.push(`Multi-sample VCF: using the first sample (${cols[9]}) of ${cols.length - 9}.`);
      }
      continue;
    }
    if (!line) continue;
    const f = line.split('\t');
    if (f.length < 10) continue;
    const fmt = f[8].split(':');
    const gtIdx = fmt.indexOf('GT');
    if (gtIdx < 0) continue;
    const gt = f[9].split(':')[gtIdx];
    if (!gt || gt === '.' || gt === './.' || gt === '.|.') continue;
    const alleles = [f[3], ...f[4].split(',')];
    let ok = true;
    const bases = [];
    for (const idx of gt.split(/[/|]/)) {
      if (idx === '.') { ok = false; break; }
      const a = alleles[parseInt(idx, 10)];
      if (a === undefined || a.length !== 1) { ok = false; break; }
      bases.push(a);
    }
    if (!ok || !bases.length) continue;
    const chrom = normalizeChromosome(f[0]);
    const pos = parseInt(f[1], 10);
    calls.push({
      rsid: f[2] && f[2] !== '.' ? f[2] : '',
      chrom,
      pos: Number.isNaN(pos) ? null : pos,
      g: normalizePloidy(chrom, normalizeGenotype(bases.join(''))),
      build,
    });
  }
  if (build === null) warnings.push('VCF declares no reference assembly.');
  return { calls, format: 'vcf', build, warnings, orientation: {}, flipped: {}, curation: {} };
}

/* Promethease reports embed base64 + zlib JSON payloads. */
async function inflate(bytes) {
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate'));
  return new Response(stream).text();
}

async function parsePromethease(text, onProgress) {
  const re = /decompressString\('([A-Za-z0-9+/=]+)'\)/g;
  const calls = [];
  const orientation = {}, flipped = {}, curation = {}, table = {};
  const builds = new Set();
  let match, n = 0;
  const CLNSIG = {1:'Untested',2:'Non-pathogenic',3:'Probable non-pathogenic',
    4:'Probable pathogenic',5:'Pathogenic',6:'Drug response',
    7:'Histocompatibility',255:'Other'};

  const payloads = [];
  while ((match = re.exec(text)) !== null) payloads.push(match[1]);

  for (const encoded of payloads) {
    if (++n % 500 === 0 && onProgress) {
      onProgress(`decoding report ${Math.round(100 * n / payloads.length)}%`);
      await new Promise(r => setTimeout(r, 0));
    }
    let chunk;
    try {
      const bin = atob(encoded);
      const bytes = new Uint8Array(bin.length);
      for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
      chunk = JSON.parse(await inflate(bytes));
    } catch (e) { continue; }

    if (Array.isArray(chunk)) {
      for (const rec of chunk) {
        if (!rec || typeof rec !== 'object' || !rec.rsnum) continue;
        if (typeof rec.reference === 'number') builds.add(rec.reference);
        const chrom = normalizeChromosome(rec.chrom);
        // `was` is the plus-strand chip call; `geno` is SNPedia-oriented.
        let plus = normalizeGenotype(rec.was);
        const oriented = normalizeGenotype(rec.geno);
        if (!plus) { plus = oriented; if (rec.flipped && plus) plus = complement(plus); }
        calls.push({
          rsid: rec.rsnum, chrom,
          pos: typeof rec.pos === 'number' ? rec.pos : null,
          g: normalizePloidy(chrom, plus),
          build: typeof rec.reference === 'number' ? rec.reference : null,
        });
        if (rec.orientation) orientation[rec.rsnum] = rec.orientation;
        flipped[rec.rsnum] = !!rec.flipped;
        curation[rec.rsnum] = {
          summary: rec.genosummary || '', repute: rec.repute || null,
          magnitude: typeof rec.magnitude === 'number' ? rec.magnitude : null,
          genes: rec.genes || [],
          clinvar: CLNSIG[rec.clinvar_1] || null,
        };
      }
    } else if (chunk && typeof chunk === 'object') {
      // rsid -> genotype -> {repute, mag}
      for (const [rsid, genos] of Object.entries(chunk)) {
        if (!genos || typeof genos !== 'object') continue;
        const bucket = table[rsid] || (table[rsid] = {});
        for (const [rawG, grade] of Object.entries(genos)) {
          if (!grade || typeof grade !== 'object') continue;
          const key = normalizeGenotype(rawG);
          if (!key) continue;
          const mag = parseFloat(grade.mag);
          bucket[key] = { repute: grade.repute || null, magnitude: Number.isNaN(mag) ? null : mag };
        }
      }
    }
  }

  const warnings = [];
  const minus = Object.values(orientation).filter(v => v === 'minus').length;
  if (minus) {
    warnings.push(`${minus.toLocaleString()} of ${Object.keys(orientation).length.toLocaleString()} ` +
      'SNPs are reported on the minus strand; orientation is preserved per SNP.');
  }
  return {
    calls, format: 'promethease',
    build: builds.size === 1 ? [...builds][0] : null,
    warnings, orientation, flipped, curation, snpediaTable: table,
  };
}

async function parseFile(text, onProgress) {
  const head = text.slice(0, 4096);
  if (head.toLowerCase().includes('<html') || head.toLowerCase().includes('promethease')) {
    return parsePromethease(text, onProgress);
  }
  if (head.trimStart().startsWith('##fileformat=VCF')) return parseVcf(text);
  return parseChipText(text);
}

/* ----------------------------------------------------------- annotation -- */

const ZYG = { ABSENT:'absent', HET:'heterozygous', HOM:'homozygous', HEMI:'hemizygous', UNKNOWN:'unknown' };
const CHIP_ERROR_RATE = 1e-3;
const COHORT_DETECTION_LIMIT = 1e-5;
const COMMON_CEILING = 0.05;
const MIN_MAGNITUDE = 2.0;

function zygosityFor(genotype, allele) {
  if (!genotype || !allele || allele.length !== 1) return ZYG.UNKNOWN;
  const count = genotype.split('').filter(c => c === allele).length;
  if (count === 0) return ZYG.ABSENT;
  if (genotype.length === 1) return ZYG.HEMI;
  return count === 2 ? ZYG.HOM : ZYG.HET;
}

function plausibility(zygosity, frequency, significance) {
  const path = significance && significance.toLowerCase().includes('pathogenic');

  if (frequency != null && frequency >= 0 && path && frequency > COMMON_CEILING) {
    return `classified pathogenic yet carried by ${(frequency * 100).toFixed(1)}% of people. ` +
      'A variant that common cannot cause a rare dominant condition, so this is a ' +
      'reference-allele or annotation artifact rather than a finding.';
  }

  if (frequency != null && frequency > 0) {
    let expected = null;
    if (zygosity === ZYG.HOM) expected = frequency * frequency;
    else if (zygosity === ZYG.HET) expected = 2 * frequency * (1 - frequency);
    else if (zygosity === ZYG.HEMI) expected = frequency;
    if (expected !== null && expected < CHIP_ERROR_RATE) {
      const oneIn = Math.round(1 / expected).toLocaleString();
      const ratio = Math.round(CHIP_ERROR_RATE / expected).toLocaleString();
      return `${zygosity} for an allele found in ${(frequency * 100).toPrecision(3)}% of people, ` +
        `expected in roughly 1 in ${oneIn}. A genotyping error is on the order of ${ratio}x ` +
        'more likely than a true call.';
    }
    return null;
  }

  // No frequency at all: absent from 1000 Genomes, ExAC and ESP because it is
  // rarer than those cohorts can resolve, not because nobody looked.
  if (path && (zygosity === ZYG.HOM || zygosity === ZYG.HEMI)) {
    const dosage = zygosity === ZYG.HOM ? 'homozygous' : 'hemizygous';
    return `${dosage} for a variant classified pathogenic that appears in no population ` +
      'frequency cohort (1000 Genomes, ExAC, ESP), placing it below roughly 1 in 100,000 ' +
      'carriers. On a consumer array a genotyping error is overwhelmingly more likely. ' +
      'Treat as an artifact unless confirmed by clinical-grade sequencing.';
  }
  return null;
}

class Bundle {
  constructor(clinvar, gwas) {
    this.clinvar = clinvar;
    this.gwas = gwas;
    this.clinvarIndex = Bundle.index(clinvar.rows);
    this.gwasIndex = Bundle.index(gwas.rows);
  }
  static index(rows) {
    const map = new Map();
    for (let i = 0; i < rows.length; i++) {
      const rs = rows[i][0];
      const at = map.get(rs);
      if (at === undefined) map.set(rs, [i]); else at.push(i);
    }
    return map;
  }
  static rsNum(rsid) {
    if (!rsid || !rsid.toLowerCase().startsWith('rs')) return null;
    const n = parseInt(rsid.slice(2), 10);
    return Number.isNaN(n) ? null : n;
  }
}

function annotate(sample, bundle, onProgress) {
  const findings = [];
  let considered = 0;

  for (const call of sample.calls) {
    if (!usable(call)) continue;
    considered++;
    const rs = Bundle.rsNum(call.rsid);
    const annotations = [];

    if (rs !== null) {
      for (const i of bundle.clinvarIndex.get(rs) || []) {
        const [, alt, sigIdx, stars, freq, condIdx, geneIdx] = bundle.clinvar.rows[i];
        const zygosity = zygosityFor(call.g, alt);
        if (zygosity === ZYG.ABSENT || zygosity === ZYG.UNKNOWN) continue;
        const significance = bundle.clinvar.significance[sigIdx] || '';
        const frequency = freq < 0 ? null : freq;
        annotations.push({
          source: 'clinvar', category: 'clinical',
          title: significance, significance,
          stars: stars < 0 ? null : stars,
          conditions: condIdx >= 0 ? [bundle.clinvar.conditions[condIdx]] : [],
          genes: geneIdx >= 0 ? [bundle.clinvar.genes[geneIdx]] : [],
          zygosity, frequency,
          flag: plausibility(zygosity, frequency, significance),
        });
      }
      for (const i of bundle.gwasIndex.get(rs) || []) {
        const [, traitIdx, risk, geneIdx] = bundle.gwas.rows[i];
        // No reported risk allele means carriage cannot be established, so
        // the association says nothing about this person. Mirrors `applies`.
        if (!risk) continue;
        const zygosity = zygosityFor(call.g, risk);
        if (zygosity === ZYG.ABSENT || zygosity === ZYG.UNKNOWN) continue;
        annotations.push({
          source: 'gwas', category: 'trait',
          title: bundle.gwas.traits[traitIdx] || '',
          genes: geneIdx >= 0 ? [bundle.gwas.genes[geneIdx]] : [],
          conditions: [], zygosity, riskAllele: risk || null,
          note: 'association only, not a diagnosis',
        });
      }
    }

    const curated = sample.curation && sample.curation[call.rsid];
    if (curated) {
      let g = call.g, ambiguous = false;
      if (sample.flipped[call.rsid]) {
        if (isStrandAmbiguous(g)) ambiguous = true; else g = complement(g);
      }
      const graded = sample.snpediaTable && sample.snpediaTable[call.rsid]
        ? sample.snpediaTable[call.rsid][g] : null;
      const magnitude = graded && graded.magnitude != null ? graded.magnitude : curated.magnitude;
      if (curated.summary || magnitude != null) {
        annotations.push({
          source: 'snpedia', category: 'curated',
          title: curated.summary || '',
          magnitude, repute: (graded && graded.repute) || curated.repute,
          genes: curated.genes || [], conditions: [],
          zygosity: call.g.length === 1 ? ZYG.HEMI
            : (call.g[0] !== call.g[1] ? ZYG.HET : ZYG.HOM),
          flag: ambiguous
            ? 'palindromic A/T or C/G SNP on a flipped entry: orientation cannot be verified from the alleles'
            : null,
        });
      }
    }

    if (!annotations.length) continue;

    // Materiality: a ClinVar entry with assertion criteria, a SNPedia magnitude
    // of 2+, or a genome-wide significant association actually carried.
    const material = annotations.some(a =>
      (a.source === 'clinvar' && (a.stars === null || a.stars >= 1)) ||
      (a.source === 'snpedia' && a.magnitude >= MIN_MAGNITUDE) ||
      (a.source === 'gwas'));
    if (!material) continue;

    const flags = annotations.map(a => a.flag).filter(Boolean);
    const conflicts = [];
    const path = annotations.some(a => (a.significance || '').toLowerCase().includes('pathogenic'));
    if (path && annotations.some(a => a.source === 'snpedia' && a.repute === 'Good')) {
      conflicts.push('SNPedia grades this Good while ClinVar classifies it pathogenic; ' +
        'the report curation is a snapshot and ClinVar is current');
    }

    findings.push({
      rsid: call.rsid, genotype: call.g, chrom: call.chrom,
      annotations, flags, conflicts,
      implausible: flags.length > 0 && annotations.some(a => a.flag && a.source === 'clinvar'),
      score: scoreOf(annotations),
    });
  }

  findings.sort((a, b) => (a.implausible - b.implausible) || (b.score - a.score));
  return { findings, considered, calls: sample.calls.length };
}

const SEVERITY = { 'pathogenic':100, 'likely pathogenic':90, 'pathogenic/likely pathogenic':95,
  'drug response':60, 'risk factor':55 };

function scoreOf(annotations) {
  let best = 0;
  for (const a of annotations) {
    if (a.source === 'clinvar') {
      const w = SEVERITY[(a.significance || '').toLowerCase()] || 40;
      const stars = a.stars === null ? 1 : a.stars;
      best = Math.max(best, w * (0.4 + 0.15 * stars));
    } else if (a.source === 'snpedia' && a.magnitude) {
      best = Math.max(best, a.magnitude * 8);
    } else if (a.source === 'gwas') {
      best = Math.max(best, 10);
    }
  }
  return best;
}

window.Allele = { parseFile, annotate, Bundle, ZYG };
