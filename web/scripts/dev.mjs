import { createServer } from "node:http";
import { createReadStream, existsSync, statSync } from "node:fs";
import { extname, join, normalize } from "node:path";

const root = new URL("../public/", import.meta.url).pathname;
const manuscriptCounts = new URL(
  "../../data/derived/GSE70038_mageck_counts.tsv",
  import.meta.url,
).pathname;
const port = Number(process.env.PORT || 4173);
const host = process.env.HOST || "127.0.0.1";
const types = {
  ".css": "text/css; charset=utf-8",
  ".csv": "text/csv; charset=utf-8",
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript; charset=utf-8",
  ".json": "application/json; charset=utf-8",
  ".svg": "image/svg+xml",
  ".tsv": "text/tab-separated-values; charset=utf-8",
};

const server = createServer((request, response) => {
  const pathname = decodeURIComponent(new URL(request.url, "http://localhost").pathname);
  const relative = pathname === "/" ? "index.html" : pathname.replace(/^\/+/, "");
  const candidate = normalize(join(root, relative));
  const file = pathname === "/manuscript-gse70038-counts.tsv"
    ? manuscriptCounts
    : candidate.startsWith(root) && existsSync(candidate) &&
      statSync(candidate).isFile() ? candidate : join(root, "index.html");
  response.writeHead(200, {
    "Content-Type": types[extname(file)] || "application/octet-stream",
    "Cache-Control": "no-store",
  });
  createReadStream(file).pipe(response);
});

server.listen(port, host, () => {
  const displayHost = host === "0.0.0.0" ? "localhost" : host;
  console.log(`Local URL: http://${displayHost}:${port}`);
});
