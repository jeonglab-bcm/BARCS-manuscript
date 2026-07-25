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
    "dist/client/examples/liang-hap1/liang-hap1-full-counts.csv",
    "dist/client/examples/liang-hap1/liang-hap1-full-guide-library.csv",
    "dist/client/examples/liang-hap1/liang-hap1-metadata.csv",
    "dist/client/examples/liang-hap1/liang-hap1-real-fastq-manifest.csv",
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
  assert.match(html, /Run full Liang HAP1/);
  assert.match(html, /Prepare real Liang FASTQs/);
  assert.match(html, /Verify manuscript/);
  assert.match(html, /app\.js\?v=20260725-stackfix/);

  const app = readFileSync("dist/client/app.js", "utf8");
  const csv = readFileSync("dist/client/csv.js", "utf8");
  const fastqWorker = readFileSync("dist/client/fastq-worker.js", "utf8");
  assert.doesNotMatch(app, /Math\.(?:max|min)\(\.\.\./);
  assert.doesNotMatch(csv, /\.\.\.rows\.map/);
  assert.doesNotMatch(fastqWorker, /push\(\.\.\.parts\)/);
});
