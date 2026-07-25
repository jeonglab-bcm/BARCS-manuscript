import {
  calibrateControls,
  geneConsistency,
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
      geneDiagnostics = geneConsistency(guideResults, data.options.gene);
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
            center: geneDiagnostics.center,
            scale: geneDiagnostics.scale,
            controlGenes: geneDiagnostics.controlGenes,
            usedControlNull: geneDiagnostics.usedControlNull,
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
