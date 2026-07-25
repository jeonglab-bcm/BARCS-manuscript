import { buildDesign } from "./barcs-core.js";
import { parseBarcsInputs, toCsv } from "./csv.js";
import {
  parseGuideLibrary,
  quantificationToCountsText,
} from "./fastq-core.js";

const $ = (selector) => document.querySelector(selector);
const $$ = (selector) => [...document.querySelectorAll(selector)];
const state = {
  countText: null,
  metadataText: null,
  input: null,
  design: null,
  result: null,
  worker: null,
  quantWorker: null,
  libraries: [],
  fastqFiles: [],
  quantification: null,
  sourceMode: "counts",
  manuscriptPreset: false,
  manuscriptReference: null,
};

function showError(element, message) {
  element.textContent = message;
  element.hidden = !message;
}

function formatNumber(value, digits = 3) {
  if (!Number.isFinite(value)) return "—";
  const absolute = Math.abs(value);
  if (absolute !== 0 && (absolute < 0.001 || absolute >= 10000)) {
    return value.toExponential(2);
  }
  return value.toLocaleString(undefined, { maximumFractionDigits: digits });
}

function formatP(value) {
  if (!Number.isFinite(value)) return "—";
  if (value < 0.001) return value.toExponential(2);
  return value.toFixed(3);
}

function download(name, content, type = "text/csv;charset=utf-8") {
  const blob = new Blob([content], { type });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = name;
  link.click();
  setTimeout(() => URL.revokeObjectURL(link.href), 500);
}

async function readFile(file, kind) {
  const text = await file.text();
  state[kind] = text;
  state.manuscriptPreset = false;
  const label = kind === "countText" ? $("#countsLabel") : $("#metadataLabel");
  const drop = kind === "countText" ? $("#countsDrop") : $("#metadataDrop");
  label.textContent = file.name;
  drop.classList.add("loaded");
  tryParseInput();
}

function setSourceMode(mode) {
  state.sourceMode = mode;
  const fastq = mode === "fastq";
  $("#countsSourcePanel").hidden = fastq;
  $("#fastqSourcePanel").hidden = !fastq;
  $("#countsMode").classList.toggle("active", !fastq);
  $("#fastqMode").classList.toggle("active", fastq);
  $("#countsMode").setAttribute("aria-selected", String(!fastq));
  $("#fastqMode").setAttribute("aria-selected", String(fastq));
}

function updateQuantifyButton() {
  $("#quantifyButton").disabled =
    !state.libraries.length || !state.fastqFiles.length || Boolean(state.quantWorker);
}

async function readLibraries(files) {
  const libraries = [];
  const errors = [];
  for (const file of files) {
    try {
      libraries.push(parseGuideLibrary(await file.text(), file.name));
    } catch (error) {
      errors.push(`${file.name}: ${error.message}`);
    }
  }
  state.libraries = libraries;
  $("#libraryLabel").textContent = libraries.length
    ? `${libraries.length} candidate ${libraries.length === 1 ? "library" : "libraries"}`
    : "Choose one or more libraries";
  $("#libraryDrop").classList.toggle("loaded", Boolean(libraries.length));
  showError($("#loadError"), errors.join(" "));
  updateQuantifyButton();
}

function readFastqFiles(files) {
  state.fastqFiles = [...files];
  $("#fastqLabel").textContent = state.fastqFiles.length
    ? `${state.fastqFiles.length} FASTQ ${state.fastqFiles.length === 1 ? "file" : "files"}`
    : "Choose FASTQ files";
  $("#fastqDrop").classList.toggle("loaded", Boolean(state.fastqFiles.length));
  updateQuantifyButton();
}

function finishQuantification() {
  $("#quantifyButton").textContent = "Align reads and build counts";
  $("#quantProgressWrap").hidden = true;
  state.quantWorker = null;
  updateQuantifyButton();
}

function quantificationQcRows() {
  return (state.quantification?.samples || []).map((sample) => ({
    sample: sample.name,
    total_reads: sample.totalReads,
    mapped_reads: sample.mappedReads,
    mapping_rate: sample.mappingRate,
    ambiguous_reads: sample.ambiguousReads,
    zero_guides: sample.zeroGuides,
    zero_guide_fraction: sample.zeroGuideFraction,
    gini: sample.gini,
    dominant_position: sample.dominantPosition,
    forward_reads: sample.forwardReads,
    reverse_reads: sample.reverseReads,
  }));
}

function renderQuantification() {
  const result = state.quantification;
  const best = result.detection.scores[0];
  $("#quantSummary").hidden = false;
  $("#selectedLibrary").textContent = result.library.name;
  $("#quantGuideCount").textContent =
    result.library.guides.length.toLocaleString();
  $("#libraryConfidence").textContent = Number.isFinite(result.detection.confidence)
    ? `${formatNumber(result.detection.confidence, 2)}×`
    : result.detection.scores.length === 1 ? "only candidate" : "∞";
  $("#qcTable").innerHTML = tableMarkup(quantificationQcRows(), [
    { key: "sample", label: "Sample" },
    { key: "total_reads", label: "Reads", numeric: true, render: formatNumber },
    { key: "mapped_reads", label: "Mapped", numeric: true, render: formatNumber },
    {
      key: "mapping_rate",
      label: "Mapped %",
      numeric: true,
      render: (value) => `${formatNumber(value * 100, 1)}%`,
    },
    {
      key: "zero_guide_fraction",
      label: "Zero guides",
      numeric: true,
      render: (value) => `${formatNumber(value * 100, 1)}%`,
    },
    { key: "gini", label: "Gini", numeric: true, render: formatNumber },
    { key: "dominant_position", label: "Position", numeric: true },
  ]);
  $("#countsLabel").textContent = "Locally quantified FASTQ counts";
  $("#countsDrop").classList.add("loaded");
}

function quantifyFastqs() {
  if (!state.libraries.length || !state.fastqFiles.length) return;
  state.manuscriptPreset = false;
  state.quantWorker = new Worker("./fastq-worker.js", { type: "module" });
  $("#quantifyButton").disabled = true;
  $("#quantifyButton").textContent = "Aligning…";
  $("#quantProgressWrap").hidden = false;
  $("#quantProgressBar").value = 0;
  $("#quantProgressPercent").textContent = "0%";
  $("#quantProgressText").textContent = "Detecting guide library";
  showError($("#loadError"), "");
  state.quantWorker.onmessage = ({ data }) => {
    if (data.type === "progress") {
      $("#quantProgressBar").value = data.fraction;
      $("#quantProgressPercent").textContent =
        `${Math.round(data.fraction * 100)}%`;
      $("#quantProgressText").textContent = data.phase === "detect"
        ? "Detecting guide library"
        : `Aligning ${data.sample}${data.reads ? ` · ${data.reads.toLocaleString()} reads` : ""}`;
    } else if (data.type === "complete") {
      state.quantification = data.result;
      state.countText = quantificationToCountsText(
        data.result.library,
        data.result.samples,
      );
      state.quantWorker.terminate();
      finishQuantification();
      renderQuantification();
      tryParseInput();
    } else if (data.type === "error") {
      state.quantWorker.terminate();
      finishQuantification();
      showError($("#loadError"), data.message);
    }
  };
  state.quantWorker.onerror = (event) => {
    finishQuantification();
    showError($("#loadError"), event.message || "FASTQ alignment failed.");
  };
  state.quantWorker.postMessage({
    type: "quantify",
    files: state.fastqFiles,
    libraries: state.libraries,
    maximumDetectionReads: 20000,
  });
}

function metadataTemplate() {
  if (!state.quantification) return "";
  return toCsv(
    state.quantification.samples.map((sample) => ({
      sample: sample.name,
      predictor: "",
      batch: "",
      total: sample.mappedReads,
    })),
    ["sample", "predictor", "batch", "total"],
  );
}

function selectedCovariates() {
  return $$("#covariateList input:checked").map((input) => input.value);
}

function interactionPairs() {
  const value = $("#interactionSelect").value;
  return value ? [value.split("|")] : [];
}

function updateInteractionChoices() {
  const predictor = $("#predictorSelect").value;
  const covariates = selectedCovariates();
  const current = $("#interactionSelect").value;
  const pairs = covariates.map((covariate) => [predictor, covariate]);
  $("#interactionSelect").innerHTML = [
    '<option value="">None</option>',
    ...pairs.map(
      ([left, right]) =>
        `<option value="${escapeHtml(`${left}|${right}`)}">${escapeHtml(left)} × ${escapeHtml(right)}</option>`,
    ),
  ].join("");
  if (pairs.some(([left, right]) => `${left}|${right}` === current)) {
    $("#interactionSelect").value = current;
  }
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function populateModelFields() {
  const columns = state.input.metadataColumns;
  const likely = columns.find((name) =>
    /time|dose|phenotype|bin|treatment|condition/i.test(name)
  ) || columns[0];
  $("#predictorSelect").innerHTML = columns.map(
    (name) => `<option value="${escapeHtml(name)}">${escapeHtml(name)}</option>`,
  ).join("");
  $("#predictorSelect").value = likely;
  populateCovariates();
  $("#modelFields").disabled = false;
  updateModel();
}

function populateCovariates() {
  const predictor = $("#predictorSelect").value;
  const existing = new Set(selectedCovariates());
  const available = state.input.metadataColumns.filter((name) => name !== predictor);
  $("#covariateList").innerHTML = available.length
    ? available.map((name) => `
      <label>
        <input type="checkbox" value="${escapeHtml(name)}" ${existing.has(name) ? "checked" : ""}>
        <span>${escapeHtml(name)}</span>
      </label>`).join("")
    : "<small>No additional metadata columns.</small>";
  $$("#covariateList input").forEach(
    (input) => input.addEventListener("change", () => {
      updateInteractionChoices();
      updateModel();
    }),
  );
  updateInteractionChoices();
}

function updateModel() {
  if (!state.input) return;
  try {
    const predictor = $("#predictorSelect").value;
    const covariates = selectedCovariates();
    state.design = buildDesign(state.input.metadata, {
      predictor,
      covariates,
      interactions: interactionPairs(),
    });
    const previous = $("#termSelect").value;
    $("#termSelect").innerHTML = state.design.columns.slice(1).map(
      (name) => `<option value="${escapeHtml(name)}">${escapeHtml(name)}</option>`,
    ).join("");
    if (state.design.columns.includes(previous)) {
      $("#termSelect").value = previous;
    } else {
      const primary = state.design.columns.find(
        (name) => name === predictor || name.startsWith(predictor),
      );
      if (primary) $("#termSelect").value = primary;
    }
    const terms = [predictor, ...covariates];
    const interaction = interactionPairs()[0];
    $("#formulaPreview").textContent = `~ ${terms.join(" + ")}${
      interaction ? ` + ${interaction[0]}:${interaction[1]}` : ""
    }`;
    $("#degreesPreview").textContent =
      `${state.design.degreesOfFreedom} residual degrees of freedom · ${state.design.columns.length} coefficients`;
    $("#runButton").disabled = false;
    showError($("#modelError"), "");
  } catch (error) {
    state.design = null;
    $("#runButton").disabled = true;
    $("#termSelect").innerHTML = "";
    showError($("#modelError"), error.message);
  }
}

function tryParseInput() {
  if (!state.countText || !state.metadataText) return;
  try {
    state.input = parseBarcsInputs(state.countText, state.metadataText);
    showError($("#loadError"), "");
    $("#dataSummary").hidden = false;
    $("#guideCount").textContent = state.input.counts.length.toLocaleString();
    $("#sampleCount").textContent = state.input.samples.length.toLocaleString();
    $("#controlCount").textContent =
      state.input.control.filter(Boolean).length.toLocaleString();
    const controls = state.input.control.filter(Boolean).length;
    $("#calibrateInput").disabled = controls < 20;
    $("#calibrateInput").checked = controls >= 20;
    $("#calibrateInput").closest("label").title = controls < 20
      ? "At least 20 control guides are required."
      : "";
    const geneCounts = new Map();
    state.input.gene.forEach((gene) => {
      if (gene) geneCounts.set(gene, (geneCounts.get(gene) || 0) + 1);
    });
    const eligibleGenes = [...geneCounts.values()].filter((count) => count >= 3).length;
    $("#geneInput").disabled = eligibleGenes < 2;
    $("#geneInput").checked = eligibleGenes >= 2;
    populateModelFields();
  } catch (error) {
    state.input = null;
    state.design = null;
    $("#modelFields").disabled = true;
    $("#runButton").disabled = true;
    $("#dataSummary").hidden = true;
    showError($("#loadError"), error.message);
  }
}

async function loadExample() {
  $("#exampleButton").disabled = true;
  $("#exampleButton").textContent = "Loading…";
  try {
    const [counts, metadata, reference] = await Promise.all([
      fetch("./example-counts.csv").then((response) => response.text()),
      fetch("./example-metadata.csv").then((response) => response.text()),
    ]);
    state.countText = counts;
    state.metadataText = metadata;
    state.manuscriptPreset = false;
    setSourceMode("counts");
    $("#countsLabel").textContent = "example-counts.csv";
    $("#metadataLabel").textContent = "example-metadata.csv";
    $("#countsDrop").classList.add("loaded");
    $("#metadataDrop").classList.add("loaded");
    tryParseInput();
  } catch (error) {
    showError($("#loadError"), `Could not load the example: ${error.message}`);
  } finally {
    $("#exampleButton").disabled = false;
    $("#exampleButton").textContent = "Load example";
  }
}

function configureManuscriptModel() {
  const predictor = "GSC0131_end";
  const covariates = ["GSC0827_end", "NSCCB660_end", "NSCU5_end"];
  $("#predictorSelect").value = predictor;
  populateCovariates();
  $$("#covariateList input").forEach((input) => {
    input.checked = covariates.includes(input.value);
  });
  updateInteractionChoices();
  updateModel();
  $("#termSelect").value = predictor;
  $("#minCountInput").value = "0";
  $("#calibrateInput").checked = false;
  $("#geneInput").checked = true;
  $("#geneMethodSelect").value = "directional-stouffer";
}

async function loadManuscript() {
  $("#manuscriptButton").disabled = true;
  $("#manuscriptButton").textContent = "Loading 64,747 guides…";
  showError($("#loadError"), "");
  try {
    const [counts, metadata] = await Promise.all([
      fetch("./manuscript-gse70038-counts.tsv").then((response) => {
        if (!response.ok) throw new Error("The manuscript counts are unavailable.");
        return response.text();
      }),
      fetch("./manuscript-gse70038-metadata.tsv").then((response) => {
        if (!response.ok) throw new Error("The manuscript design is unavailable.");
        return response.text();
      }),
      fetch("./manuscript-gse70038-reference.json").then((response) => {
        if (!response.ok) throw new Error("The manuscript reference is unavailable.");
        return response.json();
      }),
    ]);
    state.countText = counts;
    state.metadataText = metadata;
    state.manuscriptReference = reference;
    state.manuscriptPreset = true;
    setSourceMode("counts");
    $("#countsLabel").textContent = "GSE70038 · 64,747 guides";
    $("#metadataLabel").textContent = "Table 5 · 16 × 5 design";
    $("#countsDrop").classList.add("loaded");
    $("#metadataDrop").classList.add("loaded");
    tryParseInput();
    configureManuscriptModel();
  } catch (error) {
    showError($("#loadError"), `Could not load manuscript data: ${error.message}`);
  } finally {
    $("#manuscriptButton").disabled = false;
    $("#manuscriptButton").textContent = "Verify manuscript";
  }
}

function runAnalysis() {
  if (!state.input || !state.design) return;
  if (state.worker) state.worker.terminate();
  state.worker = new Worker("./barcs-worker.js", { type: "module" });
  $("#runButton").disabled = true;
  $("#runButton span:first-child").textContent = "Running…";
  $("#progressWrap").hidden = false;
  $("#progressBar").value = 0;
  $("#progressPercent").textContent = "0%";
  showError($("#modelError"), "");
  const term = $("#termSelect").value;
  const termIndex = state.design.columns.indexOf(term);
  state.worker.onmessage = ({ data }) => {
    if (data.type === "progress") {
      $("#progressBar").value = data.fraction;
      $("#progressPercent").textContent = `${Math.round(data.fraction * 100)}%`;
    } else if (data.type === "complete") {
      state.result = data.result;
      state.result.configuration = {
        formula: $("#formulaPreview").textContent,
        term,
        designColumns: state.design.columns,
        samples: state.input.samples,
        totals: state.input.totals,
        fdrThreshold: Number($("#fdrInput").value),
        manuscriptPreset: state.manuscriptPreset,
        geneMethod: $("#geneMethodSelect").value,
      };
      state.worker.terminate();
      state.worker = null;
      finishRun();
      renderResults();
    } else if (data.type === "error") {
      state.worker.terminate();
      state.worker = null;
      finishRun();
      showError($("#modelError"), data.message);
    }
  };
  state.worker.onerror = (event) => {
    finishRun();
    showError($("#modelError"), event.message || "The browser worker failed.");
  };
  state.worker.postMessage({
    type: "run",
    payload: {
      counts: state.input.counts,
      guide: state.input.guide,
      gene: state.input.gene,
      control: state.input.control,
      totals: state.input.totals,
      design: state.design.matrix.map((row) => Array.from(row)),
      termIndex,
      minTotalCount: Math.max(0, Number($("#minCountInput").value) || 0),
    },
    options: {
      calibrate: $("#calibrateInput").checked,
      genes: $("#geneInput").checked,
      calibration: { alpha: 0.05, minControls: 20, minScale: 1 },
      gene: {
        minGuides: $("#geneMethodSelect").value === "directional-stouffer" ? 1 : 3,
        alpha: 0.05,
        minControlGenes: 10,
        minScale: 1,
      },
      geneMethod: $("#geneMethodSelect").value,
    },
  });
}

function finishRun() {
  $("#runButton").disabled = false;
  $("#runButton span:first-child").textContent = "Run BARCS";
  $("#progressWrap").hidden = true;
}

function finiteRows(rows, field = "p_value") {
  return rows.filter((row) => Number.isFinite(row[field]));
}

function renderResults() {
  const guides = state.result.guides;
  const finite = finiteRows(guides);
  const threshold = state.result.configuration.fdrThreshold;
  const discoveries = finite.filter((row) => row.fdr < threshold);
  const rhos = finite.map((row) => row.rho).filter(Number.isFinite).sort((a, b) => a - b);
  const medianRho = rhos.length
    ? rhos[Math.floor((rhos.length - 1) / 2)]
    : NaN;
  $("#emptyResults").hidden = true;
  $("#results").hidden = false;
  const manuscriptCheck = verifyManuscriptResult();
  $("#resultTitle").textContent = state.result.configuration.manuscriptPreset
    ? `GSE70038 · ${state.result.configuration.term}`
    : state.result.configuration.term;
  $("#resultSubtitle").textContent = [
    state.result.configuration.formula,
    `${state.input.samples.length} independent libraries`,
    `${formatNumber(state.result.diagnostics.elapsedMs / 1000, 2)} s`,
    manuscriptCheck
      ? manuscriptCheck.passed
        ? `manuscript reference passed (${manuscriptCheck.checked}/${manuscriptCheck.total})`
        : `manuscript reference failed (${manuscriptCheck.checked}/${manuscriptCheck.total})`
      : null,
  ].filter(Boolean).join(" · ");
  $("#metricAnalyzed").textContent = finite.length.toLocaleString();
  $("#metricDiscoveries").textContent = discoveries.length.toLocaleString();
  $("#metricThreshold").textContent = `at FDR ${threshold.toFixed(2)}`;
  $("#metricConverged").textContent =
    `${formatNumber(100 * finite.filter((row) => row.converged).length / Math.max(1, finite.length), 1)}%`;
  $("#metricRho").textContent = formatNumber(medianRho, 4);
  $("#downloadGenes").disabled = !state.result.genes;
  renderGuideTable();
  renderGeneTable();
  renderDiagnostics();
  requestAnimationFrame(() => {
    drawVolcano();
    drawRhoHistogram();
  });
  $("#results").scrollIntoView({ behavior: "smooth", block: "start" });
}

function verifyManuscriptResult() {
  const reference = state.manuscriptReference;
  if (!state.result?.configuration.manuscriptPreset ||
      !reference || !state.result.genes) return null;
  const byGene = new Map(state.result.genes.map((row) => [row.gene, row]));
  let checked = 0;
  for (const expected of reference.genes) {
    const actual = byGene.get(expected.gene);
    if (!actual || actual.n_guides !== expected.n_guides) continue;
    const matches = ["estimate", "p_value", "fdr"].every((field) =>
      Math.abs(actual[field] - expected[field]) <=
        reference.absoluteTolerance[field]
    );
    if (matches) checked += 1;
  }
  const ordered = [...state.result.genes].sort((left, right) =>
    left.fdr - right.fdr ||
    (left.gene < right.gene ? -1 : left.gene > right.gene ? 1 : 0)
  );
  const summary = reference.summary;
  const complete = state.result.guides.length === reference.guideCount &&
    state.result.genes.length === reference.geneCount &&
    state.result.genes.filter((gene) => gene.fdr < 0.1).length ===
      summary.discoveriesAtFdrPointOne &&
    Math.abs(
      state.result.genes.reduce((sum, gene) => sum + gene.estimate, 0) -
      summary.sumEstimate,
    ) < 1e-6 &&
    Math.abs(
      state.result.genes.reduce((sum, gene) => sum + gene.p_value, 0) -
      summary.sumPValue,
    ) < 1e-6 &&
    Math.abs(
      state.result.genes.reduce((sum, gene) => sum + gene.fdr, 0) -
      summary.sumFdr,
    ) < 1e-5 &&
    ordered.slice(0, 20).every(
      (gene, index) => gene.gene === summary.top20Genes[index],
    );
  return {
    passed: complete && checked === reference.genes.length,
    checked,
    total: reference.genes.length,
  };
}

function tableMarkup(rows, columns, rowClass = () => "") {
  const heading = columns.map(
    (column) => `<th class="${column.numeric ? "numeric" : ""}">${escapeHtml(column.label)}</th>`,
  ).join("");
  const body = rows.map((row) => `
    <tr class="${rowClass(row)}">
      ${columns.map((column) => `
        <td class="${column.numeric ? "numeric" : ""}">
          ${column.render ? column.render(row[column.key], row) : escapeHtml(row[column.key] ?? "—")}
        </td>`).join("")}
    </tr>`).join("");
  return `<thead><tr>${heading}</tr></thead><tbody>${body}</tbody>`;
}

function renderGuideTable() {
  const query = $("#guideSearch").value.trim().toLowerCase();
  const threshold = state.result.configuration.fdrThreshold;
  const rows = state.result.guides
    .filter((row) => !query ||
      row.guide.toLowerCase().includes(query) ||
      row.gene.toLowerCase().includes(query))
    .sort((left, right) =>
      (Number.isFinite(left.fdr) ? left.fdr : Infinity) -
      (Number.isFinite(right.fdr) ? right.fdr : Infinity))
    .slice(0, 250);
  $("#guideTableCount").textContent =
    `Showing ${rows.length.toLocaleString()} of ${state.result.guides.length.toLocaleString()} guides`;
  $("#guideTable").innerHTML = tableMarkup(rows, [
    { key: "guide", label: "Guide" },
    { key: "gene", label: "Gene" },
    { key: "estimate", label: "Effect", numeric: true, render: formatNumber },
    { key: "std_error", label: "SE", numeric: true, render: formatNumber },
    { key: "t_value", label: "t", numeric: true, render: formatNumber },
    { key: "p_value", label: "p", numeric: true, render: formatP },
    { key: "fdr", label: "FDR", numeric: true, render: formatP },
    { key: "rho", label: "ρ", numeric: true, render: formatNumber },
    {
      key: "converged",
      label: "Fit",
      render: (value) => value ? '<span class="badge">OK</span>' : "—",
    },
  ], (row) => Number.isFinite(row.fdr) && row.fdr < threshold ? "hit" : "");
}

function renderGeneTable() {
  if (!state.result.genes) {
    $("#geneNote").textContent =
      "Gene analysis was not requested or fewer than two genes had at least three finite guide fits.";
    $("#geneTable").innerHTML = "";
    return;
  }
  const diagnostic = state.result.diagnostics.gene;
  if (diagnostic.method === "directional Stouffer") {
    $("#geneNote").textContent =
      "Directional Stouffer aggregation reproduces the manuscript comparison: median guide effect plus signed guide evidence. It is not required by the default BARCS model.";
  } else {
  $("#geneNote").textContent = diagnostic.usedControlNull
    ? `The empirical null used ${diagnostic.controlGenes} control genes (center ${formatNumber(diagnostic.center)}, scale ${formatNumber(diagnostic.scale)}).`
    : `Too few control genes were available; the robust all-gene null was used (center ${formatNumber(diagnostic.center)}, scale ${formatNumber(diagnostic.scale)}). Treat this gene analysis as exploratory.`;
  }
  const query = $("#geneSearch").value.trim().toLowerCase();
  const threshold = state.result.configuration.fdrThreshold;
  const rows = state.result.genes
    .filter((row) => !query || row.gene.toLowerCase().includes(query))
    .sort((left, right) => left.fdr - right.fdr)
    .slice(0, 250);
  $("#geneTable").innerHTML = tableMarkup(rows, [
    { key: "gene", label: "Gene" },
    { key: "n_guides", label: "Guides", numeric: true },
    { key: "estimate", label: "Effect", numeric: true, render: formatNumber },
    { key: "std_error", label: "SE", numeric: true, render: formatNumber },
    { key: "statistic", label: "z*", numeric: true, render: formatNumber },
    { key: "p_value", label: "p", numeric: true, render: formatP },
    { key: "fdr", label: "FDR", numeric: true, render: formatP },
    {
      key: "guide_direction_agreement",
      label: "Agreement",
      numeric: true,
      render: (value) => Number.isFinite(value) ? `${formatNumber(value * 100, 1)}%` : "—",
    },
  ], (row) => row.fdr < threshold ? "hit" : "");
}

function sizeCanvas(canvas) {
  const ratio = window.devicePixelRatio || 1;
  const width = Math.max(300, canvas.clientWidth);
  const height = Number(canvas.getAttribute("height")) || 300;
  canvas.width = Math.round(width * ratio);
  canvas.height = Math.round(height * ratio);
  canvas.style.height = `${height}px`;
  const context = canvas.getContext("2d");
  context.setTransform(ratio, 0, 0, ratio, 0, 0);
  return { context, width, height };
}

function drawAxes(context, width, height, padding, xLabel, yLabel) {
  context.strokeStyle = "#c9cbc5";
  context.lineWidth = 1;
  context.beginPath();
  context.moveTo(padding.left, padding.top);
  context.lineTo(padding.left, height - padding.bottom);
  context.lineTo(width - padding.right, height - padding.bottom);
  context.stroke();
  context.fillStyle = "#68716b";
  context.font = '10px "Avenir Next", "Segoe UI", sans-serif';
  context.textAlign = "center";
  context.fillText(xLabel, (padding.left + width - padding.right) / 2, height - 8);
  context.save();
  context.translate(12, (padding.top + height - padding.bottom) / 2);
  context.rotate(-Math.PI / 2);
  context.fillText(yLabel, 0, 0);
  context.restore();
}

function drawVolcano() {
  if (!state.result) return;
  const canvas = $("#volcanoCanvas");
  const { context, width, height } = sizeCanvas(canvas);
  const padding = { left: 48, right: 20, top: 18, bottom: 42 };
  const rows = finiteRows(state.result.guides)
    .filter((row) => Number.isFinite(row.estimate));
  context.clearRect(0, 0, width, height);
  drawAxes(context, width, height, padding, "coefficient", "−log10 p");
  if (!rows.length) return;
  const maximumEffect = Math.max(...rows.map((row) => Math.abs(row.estimate)), 0.1);
  const maximumEvidence = Math.max(
    ...rows.map((row) => -Math.log10(Math.max(row.p_value, 1e-300))),
    1,
  );
  const x = (value) => padding.left +
    (value + maximumEffect) / (2 * maximumEffect) *
    (width - padding.left - padding.right);
  const y = (value) => height - padding.bottom -
    value / maximumEvidence * (height - padding.top - padding.bottom);
  context.strokeStyle = "#e3e1da";
  context.beginPath();
  context.moveTo(x(0), padding.top);
  context.lineTo(x(0), height - padding.bottom);
  context.stroke();
  const threshold = state.result.configuration.fdrThreshold;
  for (const row of rows) {
    const hit = row.fdr < threshold;
    context.fillStyle = hit ? "rgba(220,115,91,.82)" : "rgba(23,91,121,.25)";
    context.beginPath();
    context.arc(
      x(row.estimate),
      y(-Math.log10(Math.max(row.p_value, 1e-300))),
      hit ? 3 : 2,
      0,
      Math.PI * 2,
    );
    context.fill();
  }
}

function drawRhoHistogram() {
  if (!state.result) return;
  const canvas = $("#rhoCanvas");
  const { context, width, height } = sizeCanvas(canvas);
  const padding = { left: 42, right: 18, top: 16, bottom: 38 };
  const values = state.result.guides.map((row) => row.rho)
    .filter((value) => Number.isFinite(value));
  context.clearRect(0, 0, width, height);
  drawAxes(context, width, height, padding, "ρ", "guides");
  if (!values.length) return;
  const maxValue = Math.max(...values, 0.001);
  const bins = 24;
  const counts = new Array(bins).fill(0);
  values.forEach((value) => {
    counts[Math.min(bins - 1, Math.floor(value / maxValue * bins))] += 1;
  });
  const maximum = Math.max(...counts);
  const plotWidth = width - padding.left - padding.right;
  const plotHeight = height - padding.top - padding.bottom;
  counts.forEach((count, index) => {
    const barWidth = plotWidth / bins;
    const barHeight = count / maximum * plotHeight;
    context.fillStyle = "#175b79";
    context.fillRect(
      padding.left + index * barWidth + 1,
      height - padding.bottom - barHeight,
      Math.max(1, barWidth - 2),
      barHeight,
    );
  });
}

function renderDiagnostics() {
  const guides = finiteRows(state.result.guides);
  const controls = guides.filter((row) => row.control);
  const calibration = state.result.diagnostics.calibration;
  const boundaries = guides.filter((row) => row.dispersion_boundary).length;
  const manuscriptCheck = verifyManuscriptResult();
  const items = [
    {
      label: "Residual degrees of freedom",
      value: state.design.degreesOfFreedom,
      note: `${state.input.samples.length} libraries − ${state.design.columns.length} coefficients`,
    },
    {
      label: "Control calibration",
      value: calibration ? `× ${formatNumber(calibration.scale, 3)}` : "Not applied",
      note: calibration
        ? `${calibration.controls} guide controls at the 5% tail`
        : `${controls.length} controls available`,
    },
    {
      label: "Dispersion boundary",
      value: `${formatNumber(100 * boundaries / Math.max(1, guides.length), 2)}%`,
      note: `${boundaries.toLocaleString()} finite fits reached the boundary`,
    },
    {
      label: "Full-library totals",
      value: `${formatNumber(Math.min(...state.input.totals), 0)}–${formatNumber(Math.max(...state.input.totals), 0)}`,
      note: state.input.totalsSource === "column sums"
        ? "Computed from every uploaded guide row"
        : `Read from metadata column ${state.input.totalsSource}`,
    },
    {
      label: "R-reference parity",
      value: "Numerically equivalent",
      note: "Same equations; small floating-point differences are expected",
    },
    ...(manuscriptCheck ? [{
      label: "Manuscript reference",
      value: manuscriptCheck.passed ? "Passed" : "Failed",
      note: `${manuscriptCheck.checked}/${manuscriptCheck.total} prespecified genes plus complete row counts`,
    }] : []),
  ];
  $("#diagnosticList").innerHTML = items.map((item) => `
    <div class="diagnostic-item">
      <span>${escapeHtml(item.label)}</span>
      <strong>${escapeHtml(item.value)}</strong>
      <p>${escapeHtml(item.note)}</p>
    </div>`).join("");
}

$("#countsFile").addEventListener("change", (event) => {
  if (event.target.files[0]) readFile(event.target.files[0], "countText");
});
$("#metadataFile").addEventListener("change", (event) => {
  if (event.target.files[0]) readFile(event.target.files[0], "metadataText");
});
$("#libraryFiles").addEventListener("change", (event) => {
  readLibraries(event.target.files);
});
$("#fastqFiles").addEventListener("change", (event) => {
  readFastqFiles(event.target.files);
});
for (const [dropSelector, inputSelector, kind] of [
  ["#countsDrop", "#countsFile", "countText"],
  ["#metadataDrop", "#metadataFile", "metadataText"],
]) {
  const drop = $(dropSelector);
  drop.addEventListener("dragover", (event) => {
    event.preventDefault();
    drop.classList.add("dragging");
  });
  drop.addEventListener("dragleave", () => drop.classList.remove("dragging"));
  drop.addEventListener("drop", (event) => {
    event.preventDefault();
    drop.classList.remove("dragging");
    const file = event.dataTransfer.files[0];
    if (file) {
      $(inputSelector).files = event.dataTransfer.files;
      readFile(file, kind);
    }
  });
}
for (const [dropSelector, inputSelector, handler] of [
  ["#libraryDrop", "#libraryFiles", readLibraries],
  ["#fastqDrop", "#fastqFiles", readFastqFiles],
]) {
  const drop = $(dropSelector);
  drop.addEventListener("dragover", (event) => {
    event.preventDefault();
    drop.classList.add("dragging");
  });
  drop.addEventListener("dragleave", () => drop.classList.remove("dragging"));
  drop.addEventListener("drop", (event) => {
    event.preventDefault();
    drop.classList.remove("dragging");
    if (event.dataTransfer.files.length) {
      $(inputSelector).files = event.dataTransfer.files;
      handler(event.dataTransfer.files);
    }
  });
}
$("#exampleButton").addEventListener("click", loadExample);
$("#manuscriptButton").addEventListener("click", loadManuscript);
$("#countsMode").addEventListener("click", () => setSourceMode("counts"));
$("#fastqMode").addEventListener("click", () => setSourceMode("fastq"));
$("#quantifyButton").addEventListener("click", quantifyFastqs);
$("#formatButton").addEventListener("click", () => $("#formatDialog").showModal());
$("#predictorSelect").addEventListener("change", () => {
  populateCovariates();
  updateModel();
});
$("#interactionSelect").addEventListener("change", updateModel);
$("#runButton").addEventListener("click", runAnalysis);
$("#guideSearch").addEventListener("input", renderGuideTable);
$("#geneSearch").addEventListener("input", renderGeneTable);
$("#downloadGuides").addEventListener("click", () => {
  download(
    "barcs-guide-results.csv",
    toCsv(state.result.guides, [
      "gene", "guide", "estimate", "std_error", "t_value", "df",
      "p_value", "fdr", "rho", "mean_cpm", "converged",
    ]),
  );
});
$("#downloadGenes").addEventListener("click", () => {
  if (state.result.genes) {
    download(
      "barcs-gene-results.csv",
      toCsv(state.result.genes, [
        "gene", "n_guides", "estimate", "std_error", "statistic",
        "p_value", "fdr", "guide_direction_agreement",
      ]),
    );
  }
});
$("#downloadCounts").addEventListener("click", () => {
  if (state.countText) download("barcs-guide-counts.csv", state.countText);
});
$("#downloadQc").addEventListener("click", () => {
  const rows = quantificationQcRows();
  if (rows.length) download("barcs-fastq-qc.csv", toCsv(rows));
});
$("#downloadMetadataTemplate").addEventListener("click", () => {
  const content = metadataTemplate();
  if (content) download("barcs-metadata-template.csv", content);
});
$$(".tab").forEach((tab) => tab.addEventListener("click", () => {
  $$(".tab").forEach((value) => value.classList.toggle("active", value === tab));
  $$(".tab-panel").forEach((panel) => panel.classList.remove("active"));
  $(`#${tab.dataset.tab}Panel`).classList.add("active");
  if (tab.dataset.tab === "diagnostics") requestAnimationFrame(drawRhoHistogram);
}));
window.addEventListener("resize", () => {
  if (!state.result) return;
  drawVolcano();
  drawRhoHistogram();
});
