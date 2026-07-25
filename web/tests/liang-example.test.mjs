import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { gunzipSync } from "node:zlib";

import { buildDesign } from "../public/barcs-core.js";
import { parseBarcsInputs, parseObjects } from "../public/csv.js";
import {
  parseGuideLibrary,
  quantifyFastqText,
} from "../public/fastq-core.js";

const root = new URL("../public/examples/liang-hap1/", import.meta.url);
const samples = [
  "HAP1_Day00_R1", "HAP1_Day00_R2",
  "HAP1_Day14_R1", "HAP1_Day14_R2",
];
const countsText = readFileSync(new URL("liang-hap1-counts.csv", root), "utf8");
const metadataText = readFileSync(
  new URL("liang-hap1-metadata.csv", root),
  "utf8",
);
const library = parseGuideLibrary(
  readFileSync(new URL("liang-hap1-guide-library.csv", root), "utf8"),
  "Liang HAP1 selected library",
);

test("Liang count example preserves the selected design and full-library totals", () => {
  const input = parseBarcsInputs(countsText, metadataText);
  assert.equal(input.guide.length, 72);
  assert.equal(input.control.filter(Boolean).length, 40);
  assert.deepEqual(input.samples, samples);
  assert.deepEqual(input.totals, [47033221, 47628084, 47672020, 47047711]);
  const design = buildDesign(input.metadata, {
    predictor: "day14",
    covariates: ["replicate"],
  });
  assert.equal(design.degreesOfFreedom, 1);
  assert.deepEqual(design.columns, ["(Intercept)", "day14", "replicateR2"]);
});

test("synthetic Liang FASTQs reproduce every expected selected-guide count", () => {
  const expected = parseObjects(readFileSync(
    new URL("liang-hap1-fastq-expected-counts.csv", root),
    "utf8",
  ));
  assert.equal(library.guides.length, 72);
  for (const sample of samples) {
    const fastq = gunzipSync(readFileSync(
      new URL(`${sample}.fastq.gz`, root),
    )).toString("utf8");
    const result = quantifyFastqText(fastq, library);
    assert.equal(result.mappedReads, result.totalReads);
    assert.equal(result.ambiguousReads, 0);
    assert.ok(result.forwardReads > result.reverseReads);
    assert.ok(result.reverseReads > 0);
    assert.deepEqual(
      result.counts,
      expected.rows.map((row) => Number(row[sample])),
    );
  }
});
