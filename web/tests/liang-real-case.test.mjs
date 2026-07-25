import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";

import { buildDesign, runScreen } from "../public/barcs-core.js";
import { parseBarcsInputs, parseObjects } from "../public/csv.js";

const root = new URL(
  "../public/examples/liang-hap1/",
  import.meta.url,
);
const countsText = readFileSync(
  new URL("liang-hap1-full-counts.csv", root),
  "utf8",
);
const metadataText = readFileSync(
  new URL("liang-hap1-metadata.csv", root),
  "utf8",
);

test("full Liang HAP1 example contains the complete deposited processed screen", () => {
  const input = parseBarcsInputs(countsText, metadataText);
  assert.equal(input.counts.length, 56_174);
  assert.equal(input.control.filter(Boolean).length, 997);
  assert.deepEqual(
    input.totals,
    [47_033_221, 47_628_084, 47_672_020, 47_047_711],
  );
  const design = buildDesign(input.metadata, {
    predictor: "day14",
    covariates: ["replicate"],
  });
  assert.deepEqual(
    design.columns,
    ["(Intercept)", "day14", "replicateR2"],
  );
  assert.equal(design.degreesOfFreedom, 1);

  const results = runScreen({
    counts: input.counts,
    guide: input.guide,
    gene: input.gene,
    control: input.control,
    totals: input.totals,
    design: design.matrix.map((row) => Array.from(row)),
    termIndex: design.columns.indexOf("day14"),
    minTotalCount: 10,
  });
  assert.equal(results.length, 56_174);
  assert.equal(
    results.filter((row) => Number.isFinite(row.p_value)).length,
    56_174,
  );
});

test("real Liang FASTQ manifest resolves four deposited HAP1 endpoint runs", () => {
  const library = parseObjects(readFileSync(
    new URL("liang-hap1-full-guide-library.csv", root),
    "utf8",
  ));
  assert.equal(library.rows.length, 56_322);

  const manifest = parseObjects(readFileSync(
    new URL("liang-hap1-real-fastq-manifest.csv", root),
    "utf8",
  ));
  assert.equal(manifest.rows.length, 4);
  assert.deepEqual(
    manifest.rows.map((row) => row.run_accession),
    ["SRR32105813", "SRR32105812", "SRR32105789", "SRR32105788"],
  );
  for (const row of manifest.rows) {
    assert.match(row.fastq_url, /^https:\/\/ftp\.sra\.ebi\.ac\.uk\//);
    assert.match(row.md5, /^[0-9a-f]{32}$/);
    assert.ok(Number(row.compressed_bytes) > 300_000_000);
    assert.equal(row.browser_filename, `${row.sample}.fastq.gz`);
  }
});
