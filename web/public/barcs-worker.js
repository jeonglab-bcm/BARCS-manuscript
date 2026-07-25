import {
  calibrateControls,
  geneConsistency,
  geneDirectionalStouffer,
  runScreen,
} from "./barcs-core.js";

self.onmessage = ({ data }) => {
  if (data.type !== "run") return;
  const started = performance.now();
  try {
    let guideResults = runScreen(
      data.payload,
      (fraction) => self.postMessage({ type: "progress", fraction }),
    );
    let calibration = null;
    if (data.options.calibrate) {
      calibration = calibrateControls(guideResults, data.options.calibration);
      guideResults = calibration.results;
    }
    let genes = null;
    let geneDiagnostics = null;
    if (data.options.genes) {
      geneDiagnostics = data.options.geneMethod === "directional-stouffer"
        ? geneDirectionalStouffer(guideResults, data.options.gene)
        : geneConsistency(guideResults, data.options.gene);
      genes = geneDiagnostics.results;
    }
    self.postMessage({
      type: "complete",
      result: {
        guides: guideResults,
        genes,
        diagnostics: {
          calibration: calibration && {
            scale: calibration.scale,
            alpha: calibration.alpha,
            controls: calibration.controls,
          },
          gene: geneDiagnostics && {
            method: geneDiagnostics.method || "shared-effect",
            center: geneDiagnostics.center ?? null,
            scale: geneDiagnostics.scale ?? null,
            controlGenes: geneDiagnostics.controlGenes ?? 0,
            usedControlNull: geneDiagnostics.usedControlNull ?? false,
          },
          elapsedMs: performance.now() - started,
        },
      },
    });
  } catch (error) {
    self.postMessage({
      type: "error",
      message: error instanceof Error ? error.message : String(error),
    });
  }
};
