import {
  copyFileSync,
  cpSync,
  mkdirSync,
  readFileSync,
  rmSync,
  writeFileSync,
} from "node:fs";
import { join } from "node:path";

const project = new URL("../", import.meta.url).pathname;
const source = join(project, "public");
const dist = join(project, "dist");
const client = join(dist, "client");

rmSync(dist, { recursive: true, force: true });
mkdirSync(client, { recursive: true });
cpSync(source, client, { recursive: true });
copyFileSync(
  join(project, "..", "data", "derived", "GSE70038_mageck_counts.tsv"),
  join(client, "manuscript-gse70038-counts.tsv"),
);
mkdirSync(join(dist, "server"), { recursive: true });
mkdirSync(join(dist, ".openai"), { recursive: true });

writeFileSync(
  join(dist, "server", "index.js"),
  `export default {
  async fetch(request, env) {
    if (env && env.ASSETS && typeof env.ASSETS.fetch === "function") {
      const response = await env.ASSETS.fetch(request);
      const url = new URL(request.url);
      if ((url.pathname === "/" || url.pathname === "/index.html") &&
          response.headers.get("content-type")?.includes("text/html")) {
        const headers = new Headers(response.headers);
        headers.set("content-type", "text/html; charset=utf-8");
        return new Response(
          (await response.text()).replaceAll("__BARCS_ORIGIN__", url.origin),
          { status: response.status, headers }
        );
      }
      return response;
    }
    return new Response("BARCS Web assets are unavailable.", { status: 503 });
  }
};
`,
);
writeFileSync(
  join(dist, ".openai", "hosting.json"),
  readFileSync(join(project, ".openai", "hosting.json")),
);

const required = [
  "index.html",
  "app.js",
  "format.js",
  "barcs-core.js",
  "barcs-worker.js",
  "fastq-core.js",
  "fastq-worker.js",
  "manuscript-gse70038-counts.tsv",
  "manuscript-gse70038-metadata.tsv",
  "styles.css",
];
for (const file of required) {
  readFileSync(join(client, file));
}
console.log("Built BARCS Web into dist/.");
