import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import {
  determineLibrary,
  parseGuideLibrary,
  quantifyFastqText,
  quantificationToCountsText,
  reverseComplement,
  stripFastqExtension,
} from "../public/fastq-core.js";
import { parseObjects } from "../public/csv.js";

const toyDirectory = new URL(
  "../../CB2/inst/extdata/toydata/",
  import.meta.url,
);
const fasta = readFileSync(new URL("small_sample.fasta", toyDirectory), "utf8");
const library = parseGuideLibrary(fasta, "CB2 toy library");

test("FASTA library parsing preserves CB2 guide annotations", () => {
  assert.equal(library.totalGuides, 25);
  assert.equal(library.guides.length, 25);
  assert.deepEqual(library.lengths, [20]);
  assert.equal(library.guides[0].guide, "POSCTRL_1");
  assert.equal(library.guides[0].gene, "POSCTRL");
  assert.equal(library.guides[5].control, true);
});

test("repeated guide sequences are excluded from exact alignment", () => {
  const duplicated = parseGuideLibrary(
    readFileSync(new URL("small_sample_dup.fasta", toyDirectory), "utf8"),
    "duplicated toy library",
  );
  assert.equal(duplicated.totalGuides, 27);
  assert.equal(duplicated.duplicateSequences, 2);
  assert.equal(duplicated.excludedGuides, 2);
  assert.equal(duplicated.guides.some((guide) => guide.guide === "RAB_6"), false);
  assert.equal(duplicated.guides.some((guide) => guide.guide === "RAB_7"), false);
});

test("browser quantification reproduces the CB2 toy FASTQ counts", () => {
  const expected = {
    Base1: 688,
    Base2: 608,
    Low1: 1730,
    Low2: 2001,
    High1: 703,
    High2: 659,
  };
  const samples = Object.entries(expected).map(([name, total]) => {
    const result = quantifyFastqText(
      readFileSync(new URL(`${name}.fastq`, toyDirectory), "utf8"),
      library,
    );
    assert.equal(result.totalReads, total);
    assert.equal(result.mappedReads, total);
    assert.equal(result.counts.reduce((sum, value) => sum + value, 0), total);
    assert.equal(result.ambiguousReads, 0);
    assert.equal(result.dominantPosition, 30);
    return { name, ...result };
  });
  const table = parseObjects(quantificationToCountsText(library, samples));
  assert.equal(table.rows.length, 25);
  assert.deepEqual(
    table.header.slice(0, 3),
    ["guide", "gene", "control"],
  );
  assert.deepEqual(table.header.slice(3), Object.keys(expected));
});

test("reverse-complement reads and candidate-library detection are supported", () => {
  const target = library.guides[0].sequence;
  const reverseRead = `GGG${reverseComplement(target)}TTT`;
  const result = quantifyFastqText(
    `@read\n${reverseRead}\n+\n${"I".repeat(reverseRead.length)}\n`,
    library,
  );
  assert.equal(result.mappedReads, 1);
  assert.equal(result.reverseReads, 1);
  assert.equal(result.counts[0], 1);

  const decoy = parseGuideLibrary(
    ">decoy_1\nAAAAAAAAAAAAAAAAAAAA\n>decoy_2\nCCCCCCCCCCCCCCCCCCCC\n",
    "decoy",
  );
  const detection = determineLibrary(
    [target, library.guides[1].sequence, target],
    [decoy, library],
  );
  assert.equal(detection.selected, "CB2 toy library");
  assert.equal(detection.scores[0].mapped, 3);
  assert.throws(
    () => determineLibrary(["NNNNNNNNNNNNNNNNNNNN"], [decoy, library]),
    /No exact guide match/,
  );
});

test("FASTQ filenames become stable sample identifiers", () => {
  assert.equal(stripFastqExtension("Day7_A.fastq.gz"), "Day7_A");
  assert.equal(stripFastqExtension("Day7_B.fq"), "Day7_B");
});
