import assert from "node:assert/strict";
import { existsSync, readFileSync } from "node:fs";
import test from "node:test";

test("site build contains the analysis runtime and Cloudflare entry", () => {
  for (const path of [
    "dist/client/index.html",
    "dist/client/app.js",
    "dist/client/format.js",
    "dist/client/barcs-core.js",
    "dist/client/barcs-worker.js",
    "dist/client/fastq-core.js",
    "dist/client/fastq-worker.js",
    "dist/client/manuscript-gse70038-counts.tsv",
    "dist/client/manuscript-gse70038-metadata.tsv",
    "dist/client/manuscript-gse70038-reference.json",
    "dist/client/examples/liang-hap1/liang-hap1-counts.csv",
    "dist/client/examples/liang-hap1/liang-hap1-metadata.csv",
    "dist/client/examples/liang-hap1/liang-hap1-fastq-example.zip",
    "dist/server/index.js",
    "dist/.openai/hosting.json",
  ]) {
    assert.equal(existsSync(path), true, `${path} is missing`);
  }
  const html = readFileSync("dist/client/index.html", "utf8");
  assert.match(html, /Run BARCS/);
  assert.match(html, /Runs on this device/);
  assert.match(html, /not promised to be bit-for-bit identical/);
  assert.match(html, /FASTQ → counts/);
  assert.match(html, /Load Liang HAP1/);
  assert.match(html, /Run Liang FASTQ demo/);
  assert.match(html, /Verify manuscript/);
});
