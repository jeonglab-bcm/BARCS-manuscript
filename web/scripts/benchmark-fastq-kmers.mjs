import { readFileSync } from "node:fs";
import { performance } from "node:perf_hooks";
import {
  addReadToQuantification,
  buildLibraryIndex,
  createQuantification,
  parseGuideLibrary,
} from "../public/fastq-core.js";

const requestedReads = Number.parseInt(process.argv[2] || "1000000", 10);
if (!Number.isSafeInteger(requestedReads) || requestedReads < 1) {
  throw new Error("Read count must be a positive integer.");
}

const libraryPath = new URL(
  "../public/examples/liang-hap1/liang-hap1-full-guide-library.csv",
  import.meta.url,
);
const library = parseGuideLibrary(
  readFileSync(libraryPath, "utf8"),
  "Liang full library",
);

const indexStarted = performance.now();
const index = buildLibraryIndex(library);
const indexMilliseconds = performance.now() - indexStarted;
const result = createQuantification(library);

const countingStarted = performance.now();
for (let readIndex = 0; readIndex < requestedReads; readIndex += 1) {
  const guide = library.guides[readIndex % library.guides.length];
  addReadToQuantification(`GATTACA${guide.sequence}TGC`, index, result);
}
const countingSeconds = (performance.now() - countingStarted) / 1000;

console.log(JSON.stringify({
  guide_count: library.guides.length,
  guide_length: library.lengths[0],
  index_milliseconds: Number(indexMilliseconds.toFixed(1)),
  reads: requestedReads,
  counting_seconds: Number(countingSeconds.toFixed(3)),
  reads_per_second: Math.round(requestedReads / countingSeconds),
  mapped_reads: result.mappedReads,
}, null, 2));
