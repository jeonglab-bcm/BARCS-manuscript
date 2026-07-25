import assert from "node:assert/strict";
import test from "node:test";

import { formatNumber, formatP } from "../public/format.js";

test("formatNumber accepts explicit precision", () => {
  assert.equal(formatNumber(1.23456, 2), "1.23");
  assert.equal(formatNumber(1.9, 0), "2");
});

test("formatNumber ignores the table row passed as its second argument", () => {
  assert.doesNotThrow(() => formatNumber(1.23456, { guide: "g1" }));
  assert.equal(formatNumber(1.23456, { guide: "g1" }), "1.235");
});

test("formatters handle tails and non-finite values", () => {
  assert.equal(formatNumber(0.00001), "1.00e-5");
  assert.equal(formatNumber(Number.NaN), "—");
  assert.equal(formatP(0.00001), "1.00e-5");
  assert.equal(formatP(Number.NaN), "—");
});
