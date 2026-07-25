import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { buildDesign } from "../public/barcs-core.js";
import { parseBarcsInputs } from "../public/csv.js";

test("browser worker completes the example screen end to end", async () => {
  const input = parseBarcsInputs(
    readFileSync(new URL("../public/example-counts.csv", import.meta.url), "utf8"),
    readFileSync(new URL("../public/example-metadata.csv", import.meta.url), "utf8"),
  );
  const design = buildDesign(input.metadata, {
    predictor: "time",
    covariates: ["batch"],
    interactions: [],
  });
  const messages = [];
  globalThis.self = {
    postMessage(message) {
      messages.push(message);
    },
  };
  await import(`../public/barcs-worker.js?worker-test=${Date.now()}`);
  globalThis.self.onmessage({
    data: {
      type: "run",
      payload: {
        counts: input.counts,
        guide: input.guide,
        gene: input.gene,
        control: input.control,
        totals: input.totals,
        design: design.matrix.map((row) => Array.from(row)),
        termIndex: 1,
        minTotalCount: 10,
      },
      options: {
        calibrate: true,
        genes: true,
        calibration: { minControls: 20 },
        gene: { minGuides: 3 },
      },
    },
  });
  const completion = messages.find((message) => message.type === "complete");
  assert.ok(completion, "worker did not send a completion message");
  assert.equal(completion.result.guides.length, 48);
  assert.equal(completion.result.genes.length, 9);
  assert.equal(completion.result.diagnostics.calibration.controls, 24);
  assert.ok(completion.result.guides.every((row) => Number.isFinite(row.fdr)));
  delete globalThis.self;
});
