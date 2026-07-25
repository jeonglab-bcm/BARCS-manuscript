import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import {
  bhAdjust,
  buildDesign,
  calibrateControls,
  fitGuide,
  geneConsistency,
  runScreen,
  studentTwoSidedP,
} from "../public/barcs-core.js";
import {
  parseBarcsInputs,
  parseDelimited,
  parseObjects,
  toCsv,
} from "../public/csv.js";

const countsText = readFileSync(
  new URL("../public/example-counts.csv", import.meta.url),
  "utf8",
);
const metadataText = readFileSync(
  new URL("../public/example-metadata.csv", import.meta.url),
  "utf8",
);

function close(actual, expected, tolerance, label) {
  assert.ok(
    Math.abs(actual - expected) <= tolerance,
    `${label}: ${actual} differs from ${expected} by more than ${tolerance}`,
  );
}

const parityTolerance = Object.freeze({
  estimate: 1e-7,
  stdError: 2e-7,
  statistic: 2e-5,
  probability: 5e-8,
  rho: 1e-8,
  geneProbability: 3e-7,
});

function referenceRows(name) {
  return parseObjects(readFileSync(
    new URL(`fixtures/${name}`, import.meta.url),
    "utf8",
  )).rows;
}

function asNumber(value) {
  return value === "" ? NaN : Number(value);
}

test("CSV parser supports commas, tabs, and quoted fields", () => {
  assert.deepEqual(parseDelimited('a,b\n"x,y",2\n'), [["a", "b"], ["x,y", "2"]]);
  assert.deepEqual(parseDelimited("a\tb\n1\t2\n"), [["a", "b"], ["1", "2"]]);
});

test("BARCS input parser preserves full-library totals", () => {
  const input = parseBarcsInputs(countsText, metadataText);
  assert.equal(input.counts.length, 48);
  assert.equal(input.samples.length, 8);
  assert.equal(input.control.filter(Boolean).length, 24);
  assert.deepEqual(
    input.totals,
    [15547, 15373, 15266, 15124, 15076, 14963, 15059, 14952],
  );
});

test("metadata totals override filtered count-table column sums", () => {
  const parsed = parseBarcsInputs(
    "guide,gene,control,S1,S2\ng1,A,false,3,5\n",
    "sample,time,total\nS1,0,100\nS2,1,120\n",
  );
  assert.deepEqual(parsed.totals, [100, 120]);
  assert.equal(parsed.totalsSource, "total");
  assert.deepEqual(parsed.metadataColumns, ["time"]);
});

test("design construction matches additive R treatment coding", () => {
  const input = parseBarcsInputs(countsText, metadataText);
  const design = buildDesign(input.metadata, {
    predictor: "time",
    covariates: ["batch"],
    interactions: [],
  });
  assert.deepEqual(design.columns, ["(Intercept)", "time", "batchB"]);
  assert.equal(design.degreesOfFreedom, 5);
  assert.deepEqual(Array.from(design.matrix[0]), [1, 0, 0]);
  assert.deepEqual(Array.from(design.matrix[1]), [1, 0, 1]);
});

test("browser fit agrees with the reference R implementation", () => {
  const input = parseBarcsInputs(countsText, metadataText);
  const design = buildDesign(input.metadata, {
    predictor: "time",
    covariates: [],
    interactions: [],
  });
  const reference = [
    {
      estimate: -0.288560234224835,
      stdError: 0.018228981456753,
      statistic: -15.8297508234031,
      pValue: 4.03158550172796e-6,
      rho: 3.97921920586975e-6,
    },
    {
      estimate: -0.300840733470599,
      stdError: 0.0199004643974444,
      statistic: -15.1172720124679,
      pValue: 5.28334537862412e-6,
      rho: 1.09073406929713e-5,
    },
    {
      estimate: -0.294641316814602,
      stdError: 0.0173035761506007,
      statistic: -17.0277701123865,
      pValue: 2.62422295851482e-6,
      rho: 0,
    },
  ];
  reference.forEach((expected, index) => {
    const fit = fitGuide(input.counts[index], input.totals, design.matrix);
    close(fit.coefficient[1], expected.estimate, 1e-8, "coefficient");
    close(fit.standardError[1], expected.stdError, 2e-8, "standard error");
    close(fit.statistic[1], expected.statistic, 2e-6, "t statistic");
    close(fit.pValue[1], expected.pValue, 3e-12, "p-value");
    close(fit.rho, expected.rho, 2e-10, "rho");
    assert.equal(fit.degreesOfFreedom, 6);
    assert.equal(fit.converged, true);
  });
});

test("all guide results match R across supported model shapes", () => {
  const input = parseBarcsInputs(countsText, metadataText);
  const specifications = {
    time: {
      predictor: "time",
      covariates: [],
      interactions: [],
    },
    time_batch: {
      predictor: "time",
      covariates: ["batch"],
      interactions: [],
    },
    time_by_batch: {
      predictor: "time",
      covariates: ["batch"],
      interactions: [["time", "batch"]],
    },
    batch: {
      predictor: "batch",
      covariates: [],
      interactions: [],
    },
  };
  const reference = referenceRows("r-guide-reference.csv");
  for (const [model, configuration] of Object.entries(specifications)) {
    const expectedRows = reference.filter((row) => row.model === model);
    const design = buildDesign(input.metadata, configuration);
    const term = expectedRows[0].term;
    const termIndex = design.columns.indexOf(term);
    assert.notEqual(termIndex, -1, `${model} did not construct R term ${term}`);
    const raw = runScreen({
      ...input,
      design: design.matrix,
      termIndex,
    });
    const calibrated = calibrateControls(raw);
    assert.equal(raw.length, expectedRows.length);
    assert.equal(calibrated.results.length, expectedRows.length);
    expectedRows.forEach((expected, index) => {
      const actualRaw = raw[index];
      const actualCalibrated = calibrated.results[index];
      assert.equal(actualRaw.guide, expected.guide);
      assert.equal(actualRaw.gene, expected.gene);
      assert.equal(actualRaw.converged, expected.converged === "TRUE");
      close(
        actualRaw.estimate,
        asNumber(expected.estimate),
        parityTolerance.estimate,
        `${model}/${expected.guide} estimate`,
      );
      close(
        actualRaw.std_error,
        asNumber(expected.raw_std_error),
        parityTolerance.stdError,
        `${model}/${expected.guide} raw SE`,
      );
      close(
        actualRaw.t_value,
        asNumber(expected.raw_t_value),
        parityTolerance.statistic,
        `${model}/${expected.guide} raw t`,
      );
      close(
        actualRaw.p_value,
        asNumber(expected.raw_p_value),
        parityTolerance.probability,
        `${model}/${expected.guide} raw p`,
      );
      close(
        actualRaw.fdr,
        asNumber(expected.raw_fdr),
        parityTolerance.probability,
        `${model}/${expected.guide} raw FDR`,
      );
      close(
        actualRaw.rho,
        asNumber(expected.rho),
        parityTolerance.rho,
        `${model}/${expected.guide} rho`,
      );
      close(
        actualCalibrated.std_error,
        asNumber(expected.calibrated_std_error),
        parityTolerance.stdError,
        `${model}/${expected.guide} calibrated SE`,
      );
      close(
        actualCalibrated.t_value,
        asNumber(expected.calibrated_t_value),
        parityTolerance.statistic,
        `${model}/${expected.guide} calibrated t`,
      );
      close(
        actualCalibrated.p_value,
        asNumber(expected.calibrated_p_value),
        parityTolerance.probability,
        `${model}/${expected.guide} calibrated p`,
      );
      close(
        actualCalibrated.fdr,
        asNumber(expected.calibrated_fdr),
        parityTolerance.probability,
        `${model}/${expected.guide} calibrated FDR`,
      );
      assert.equal(
        actualRaw.fdr < 0.1,
        asNumber(expected.raw_fdr) < 0.1,
        `${model}/${expected.guide} raw FDR decision differs`,
      );
      assert.equal(
        actualCalibrated.fdr < 0.1,
        asNumber(expected.calibrated_fdr) < 0.1,
        `${model}/${expected.guide} calibrated FDR decision differs`,
      );
    });
    close(
      calibrated.scale,
      asNumber(expectedRows[0].control_scale),
      parityTolerance.probability,
      `${model} control scale`,
    );
  }
});

test("shared-effect gene results match the R reference", () => {
  const input = parseBarcsInputs(countsText, metadataText);
  const design = buildDesign(input.metadata, {
    predictor: "time",
    covariates: ["batch"],
    interactions: [],
  });
  const raw = runScreen({
    ...input,
    design: design.matrix,
    termIndex: design.columns.indexOf("time"),
  });
  const calibrated = calibrateControls(raw);
  const actual = geneConsistency(calibrated.results).results;
  const actualByGene = new Map(actual.map((row) => [row.gene, row]));
  const reference = referenceRows("r-gene-reference.csv");
  assert.equal(actual.length, reference.length);
  assert.deepEqual(
    actual.map((row) => row.gene),
    reference.map((row) => row.gene),
    "browser gene order differs from R",
  );
  for (const expected of reference) {
    const row = actualByGene.get(expected.gene);
    assert.ok(row, `browser gene results omitted ${expected.gene}`);
    assert.equal(row.n_guides, asNumber(expected.n_guides));
    assert.equal(row.control_gene, expected.control_gene === "TRUE");
    close(
      row.estimate,
      asNumber(expected.estimate),
      parityTolerance.estimate,
      `${expected.gene} gene estimate`,
    );
    close(
      row.std_error,
      asNumber(expected.std_error),
      parityTolerance.stdError,
      `${expected.gene} gene SE`,
    );
    close(
      row.raw_statistic,
      asNumber(expected.raw_statistic),
      parityTolerance.statistic,
      `${expected.gene} raw gene statistic`,
    );
    close(
      row.guide_direction_agreement,
      asNumber(expected.guide_direction_agreement),
      parityTolerance.probability,
      `${expected.gene} direction agreement`,
    );
    close(
      row.statistic,
      asNumber(expected.statistic),
      parityTolerance.statistic,
      `${expected.gene} calibrated gene statistic`,
    );
    close(
      row.p_value,
      asNumber(expected.p_value),
      parityTolerance.geneProbability,
      `${expected.gene} gene p`,
    );
    close(
      row.fdr,
      asNumber(expected.fdr),
      parityTolerance.geneProbability,
      `${expected.gene} gene FDR`,
    );
    assert.equal(
      row.fdr < 0.1,
      asNumber(expected.fdr) < 0.1,
      `${expected.gene} gene FDR decision differs`,
    );
  }
});

test("Student t probability and BH correction are numerically sound", () => {
  close(studentTwoSidedP(2, 10), 0.0733880347707404, 1e-12, "Student t p");
  const adjusted = bhAdjust([0.01, 0.04, 0.03, NaN]);
  close(adjusted[0], 0.03, 1e-15, "BH first");
  close(adjusted[1], 0.04, 1e-15, "BH second");
  close(adjusted[2], 0.04, 1e-15, "BH third");
  assert.equal(Number.isNaN(adjusted[3]), true);
});

test("control calibration and gene statistics retain explicit assumptions", () => {
  const guides = Array.from({ length: 30 }, (_, index) => ({
    guide: `g${index}`,
    gene: index < 15 ? `gene${Math.floor(index / 3)}` : "NTC",
    control: index >= 15,
    estimate: index < 15 ? (index % 2 ? -0.2 : 0.2) : 0.001 * (index - 22),
    std_error: 0.1,
    t_value: index < 15 ? (index % 2 ? -2 : 2) : 0.01 * (index - 22),
    df: 8,
    p_value: 0.1,
    fdr: 0.2,
    converged: true,
  }));
  const calibration = calibrateControls(guides, { minControls: 10 });
  assert.equal(calibration.scale, 1);
  const genes = geneConsistency(calibration.results, {
    minGuides: 3,
    minControlGenes: 2,
  });
  assert.ok(genes.results.length >= 6);
  assert.ok(genes.results.every((row) => "estimate" in row && "statistic" in row));
});

test("CSV export quotes identifiers without changing numeric values", () => {
  const output = toCsv([{ gene: "A,B", estimate: 1.25 }], ["gene", "estimate"]);
  assert.equal(output, 'gene,estimate\n"A,B",1.25');
});
