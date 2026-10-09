// Evaluate the browser implementation's rules over a set of cases supplied on
// stdin, and print the verdicts as JSON. `tests/test_parity.py` feeds the same
// cases through the Python engine and asserts the two agree.
//
// The comparison is rule-level rather than report-level on purpose: the CLI
// annotates against the full sqlite database and the browser against a shipped
// subset, so identical findings are not expected. Identical *rules* are.
import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';

const REPO = path.resolve(import.meta.dirname, '..', '..');

const sandbox = {
  console,
  window: {},
  document: {
    createElement: () => ({
      set innerHTML(v) { this._v = v; },
      get value() {
        return (this._v || '').replace(/&amp;/g, '&').replace(/&lt;/g, '<')
          .replace(/&gt;/g, '>').replace(/&quot;/g, '"').replace(/&#39;/g, "'");
      },
    }),
  },
};
vm.createContext(sandbox);
vm.runInContext(fs.readFileSync(path.join(REPO, 'docs', 'allele.js'), 'utf8'), sandbox);

// allele.js keeps its rule helpers module-private, so reach them through the
// same context rather than exporting them only for the tests.
const evalIn = (expr) => vm.runInContext(expr, sandbox);

const cases = JSON.parse(fs.readFileSync(0, 'utf8'));
const out = {};

out.plausibility = cases.plausibility.map((c) =>
  evalIn(`plausibility(${JSON.stringify(c.zygosity)}, ${JSON.stringify(c.frequency)},
    ${JSON.stringify(c.significance)})`) !== null);

out.zygosity = cases.zygosity.map((c) =>
  evalIn(`zygosityFor(${JSON.stringify(c.genotype)}, ${JSON.stringify(c.allele)})`));

out.normalizeGenotype = cases.normalizeGenotype.map((raw) =>
  evalIn(`normalizeGenotype(${JSON.stringify(raw)})`));

out.normalizeChromosome = cases.normalizeChromosome.map((raw) =>
  evalIn(`normalizeChromosome(${JSON.stringify(raw)})`) || null);

out.normalizePloidy = cases.normalizePloidy.map((c) =>
  evalIn(`normalizePloidy(${JSON.stringify(c.chrom)}, ${JSON.stringify(c.genotype)})`));

out.score = cases.score.map((annotations) =>
  evalIn(`scoreOf(${JSON.stringify(annotations)})`));

out.strandAmbiguous = cases.strandAmbiguous.map((g) =>
  evalIn(`isStrandAmbiguous(${JSON.stringify(g)})`));

process.stdout.write(JSON.stringify(out));
