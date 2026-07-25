import assert from "node:assert/strict";
import { existsSync, readFileSync } from "node:fs";
import test from "node:test";

test("site build contains the analysis runtime and Cloudflare entry", () => {
  for (const path of [
    "dist/client/index.html",
    "dist/client/app.js",
    "dist/client/barcs-core.js",
    "dist/client/barcs-worker.js",
    "dist/server/index.js",
    "dist/.openai/hosting.json",
  ]) {
    assert.equal(existsSync(path), true, `${path} is missing`);
  }
  const html = readFileSync("dist/client/index.html", "utf8");
  assert.match(html, /Run BARCS/);
  assert.match(html, /Runs on this device/);
});
