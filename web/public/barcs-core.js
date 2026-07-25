const EPSILON = Number.EPSILON;

export function clamp(value, lower, upper) {
  return Math.min(Math.max(value, lower), upper);
}

export function logit(value) {
  return Math.log(value / (1 - value));
}

export function logistic(value) {
  if (value >= 0) {
    const z = Math.exp(-value);
    return 1 / (1 + z);
  }
  const z = Math.exp(value);
  return z / (1 + z);
}

function logGamma(value) {
  const coefficients = [
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.6150291621406,
    12.507343278686905,
    -0.13857109526572012,
    9.984369578019572e-6,
    1.5056327351493116e-7,
  ];
  if (value < 0.5) {
    return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * value)) -
      logGamma(1 - value);
  }
  const shifted = value - 1;
  let series = 0.9999999999998099;
  for (let index = 0; index < coefficients.length; index += 1) {
    series += coefficients[index] / (shifted + index + 1);
  }
  const t = shifted + coefficients.length - 0.5;
  return 0.5 * Math.log(2 * Math.PI) +
    (shifted + 0.5) * Math.log(t) - t + Math.log(series);
}

function betaContinuedFraction(a, b, x) {
  const maxIterations = 200;
  const minimum = 1e-300;
  let qab = a + b;
  let qap = a + 1;
  let qam = a - 1;
  let c = 1;
  let d = 1 - (qab * x) / qap;
  if (Math.abs(d) < minimum) d = minimum;
  d = 1 / d;
  let result = d;
  for (let iteration = 1; iteration <= maxIterations; iteration += 1) {
    const twice = 2 * iteration;
    let numerator = iteration * (b - iteration) * x /
      ((qam + twice) * (a + twice));
    d = 1 + numerator * d;
    if (Math.abs(d) < minimum) d = minimum;
    c = 1 + numerator / c;
    if (Math.abs(c) < minimum) c = minimum;
    d = 1 / d;
    result *= d * c;

    numerator = -(a + iteration) * (qab + iteration) * x /
      ((a + twice) * (qap + twice));
    d = 1 + numerator * d;
    if (Math.abs(d) < minimum) d = minimum;
    c = 1 + numerator / c;
    if (Math.abs(c) < minimum) c = minimum;
    d = 1 / d;
    const delta = d * c;
    result *= delta;
    if (Math.abs(delta - 1) < 3e-14) break;
  }
  return result;
}

export function regularizedBeta(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const front = Math.exp(
    logGamma(a + b) - logGamma(a) - logGamma(b) +
    a * Math.log(x) + b * Math.log1p(-x),
  );
  if (x < (a + 1) / (a + b + 2)) {
    return front * betaContinuedFraction(a, b, x) / a;
  }
  return 1 - front * betaContinuedFraction(b, a, 1 - x) / b;
}

export function studentTwoSidedP(statistic, degreesOfFreedom) {
  if (!Number.isFinite(statistic) || !(degreesOfFreedom > 0)) return NaN;
  const x = degreesOfFreedom /
    (degreesOfFreedom + statistic * statistic);
  return clamp(
    regularizedBeta(x, degreesOfFreedom / 2, 0.5),
    0,
    1,
  );
}

function erf(value) {
  const sign = value < 0 ? -1 : 1;
  const x = Math.abs(value);
  const t = 1 / (1 + 0.3275911 * x);
  const polynomial = (((((1.061405429 * t - 1.453152027) * t) +
    1.421413741) * t - 0.284496736) * t + 0.254829592) * t;
  return sign * (1 - polynomial * Math.exp(-x * x));
}

export function normalCdf(value) {
  return 0.5 * (1 + erf(value / Math.SQRT2));
}

export function normalTwoSidedP(statistic) {
  return clamp(2 * (1 - normalCdf(Math.abs(statistic))), 0, 1);
}

export function normalQuantile(probability) {
  if (probability <= 0) return -Infinity;
  if (probability >= 1) return Infinity;
  const a = [
    -39.69683028665376, 220.9460984245205, -275.9285104469687,
    138.357751867269, -30.66479806614716, 2.506628277459239,
  ];
  const b = [
    -54.47609879822406, 161.5858368580409, -155.6989798598866,
    66.80131188771972, -13.28068155288572,
  ];
  const c = [
    -0.007784894002430293, -0.3223964580411365,
    -2.400758277161838, -2.549732539343734,
    4.374664141464968, 2.938163982698783,
  ];
  const d = [
    0.007784695709041462, 0.3224671290700398,
    2.445134137142996, 3.754408661907416,
  ];
  const lower = 0.02425;
  const upper = 1 - lower;
  if (probability < lower) {
    const q = Math.sqrt(-2 * Math.log(probability));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q +
      c[4]) * q + c[5]) /
      ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
  if (probability > upper) {
    const q = Math.sqrt(-2 * Math.log(1 - probability));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q +
      c[4]) * q + c[5]) /
      ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
  const q = probability - 0.5;
  const r = q * q;
  return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r +
    a[4]) * r + a[5]) * q /
    (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r +
      b[4]) * r + 1);
}

export function quantileType8(values, probability) {
  const sorted = values.filter(Number.isFinite).slice().sort((a, b) => a - b);
  if (!sorted.length) return NaN;
  if (probability <= 0) return sorted[0];
  if (probability >= 1) return sorted[sorted.length - 1];
  const h = sorted.length * probability + (probability + 1) / 3;
  const j = Math.floor(h);
  const gamma = h - j;
  const lower = sorted[clamp(j - 1, 0, sorted.length - 1)];
  const upper = sorted[clamp(j, 0, sorted.length - 1)];
  return (1 - gamma) * lower + gamma * upper;
}

export function median(values) {
  return quantileType8(values, 0.5);
}

export function bhAdjust(values) {
  const result = new Array(values.length).fill(NaN);
  const finite = values
    .map((value, index) => ({ value, index }))
    .filter(({ value }) => Number.isFinite(value))
    .sort((left, right) => left.value - right.value);
  let running = 1;
  for (let index = finite.length - 1; index >= 0; index -= 1) {
    const adjusted = finite[index].value * finite.length / (index + 1);
    running = Math.min(running, adjusted);
    result[finite[index].index] = clamp(running, 0, 1);
  }
  return result;
}

function zeros(rows, columns) {
  return Array.from({ length: rows }, () => new Float64Array(columns));
}

function weightedSystem(design, weight, response) {
  const rows = design.length;
  const columns = design[0].length;
  const information = zeros(columns, columns);
  const target = new Float64Array(columns);
  for (let row = 0; row < rows; row += 1) {
    const w = weight[row];
    const x = design[row];
    for (let left = 0; left < columns; left += 1) {
      target[left] += x[left] * w * response[row];
      for (let right = 0; right <= left; right += 1) {
        information[left][right] += x[left] * w * x[right];
      }
    }
  }
  for (let left = 0; left < columns; left += 1) {
    for (let right = 0; right < left; right += 1) {
      information[right][left] = information[left][right];
    }
  }
  return { information, target };
}

function cholesky(matrix) {
  const size = matrix.length;
  const lower = zeros(size, size);
  for (let row = 0; row < size; row += 1) {
    for (let column = 0; column <= row; column += 1) {
      let value = matrix[row][column];
      for (let index = 0; index < column; index += 1) {
        value -= lower[row][index] * lower[column][index];
      }
      if (row === column) {
        if (!(value > 1e-14) || !Number.isFinite(value)) {
          throw new Error("The design or weighted fit is singular.");
        }
        lower[row][column] = Math.sqrt(value);
      } else {
        lower[row][column] = value / lower[column][column];
      }
    }
  }
  return lower;
}

function solveFromCholesky(lower, target) {
  const size = lower.length;
  const forward = new Float64Array(size);
  const solution = new Float64Array(size);
  for (let row = 0; row < size; row += 1) {
    let value = target[row];
    for (let index = 0; index < row; index += 1) {
      value -= lower[row][index] * forward[index];
    }
    forward[row] = value / lower[row][row];
  }
  for (let row = size - 1; row >= 0; row -= 1) {
    let value = forward[row];
    for (let index = row + 1; index < size; index += 1) {
      value -= lower[index][row] * solution[index];
    }
    solution[row] = value / lower[row][row];
  }
  return solution;
}

function solvePositiveDefinite(matrix, target) {
  return solveFromCholesky(cholesky(matrix), target);
}

function inversePositiveDefinite(matrix) {
  const lower = cholesky(matrix);
  const size = matrix.length;
  const inverse = zeros(size, size);
  for (let column = 0; column < size; column += 1) {
    const unit = new Float64Array(size);
    unit[column] = 1;
    const solution = solveFromCholesky(lower, unit);
    for (let row = 0; row < size; row += 1) {
      inverse[row][column] = solution[row];
    }
  }
  return inverse;
}

function linearPredictor(design, coefficient) {
  const result = new Float64Array(design.length);
  for (let row = 0; row < design.length; row += 1) {
    let value = 0;
    for (let column = 0; column < coefficient.length; column += 1) {
      value += design[row][column] * coefficient[column];
    }
    result[row] = value;
  }
  return result;
}

function binomialInitial(count, total, design, tolerance) {
  const columns = design[0].length;
  const coefficient = new Float64Array(columns);
  const pooled = (count.reduce((sum, value) => sum + value, 0) + 0.5) /
    (total.reduce((sum, value) => sum + value, 0) + 1);
  coefficient[0] = logit(pooled);
  let beta = coefficient;
  for (let iteration = 0; iteration < 50; iteration += 1) {
    const eta = linearPredictor(design, beta);
    const response = new Float64Array(count.length);
    const weight = new Float64Array(count.length);
    for (let index = 0; index < count.length; index += 1) {
      const mu = clamp(logistic(eta[index]), 1e-8, 1 - 1e-8);
      response[index] = eta[index] +
        (count[index] / total[index] - mu) / (mu * (1 - mu));
      weight[index] = total[index] * mu * (1 - mu);
    }
    const system = weightedSystem(design, weight, response);
    const next = solvePositiveDefinite(system.information, system.target);
    let change = 0;
    for (let index = 0; index < beta.length; index += 1) {
      change = Math.max(
        change,
        Math.abs(next[index] - beta[index]) / Math.max(1, Math.abs(beta[index])),
      );
    }
    beta = next;
    if (change < tolerance) break;
  }
  return beta;
}

function estimateRho(count, total, mu, degreesOfFreedom) {
  const pearson = (rho) => {
    let result = 0;
    for (let index = 0; index < count.length; index += 1) {
      const variance = total[index] * mu[index] * (1 - mu[index]);
      const residual = count[index] - total[index] * mu[index];
      result += residual * residual /
        (variance * (1 + (total[index] - 1) * rho));
    }
    return result;
  };
  const q0 = pearson(0);
  if (!Number.isFinite(q0)) {
    return { rho: 0, scale: 1, pearson: q0, boundary: true };
  }
  if (q0 <= degreesOfFreedom) {
    return { rho: 0, scale: 1, pearson: q0, boundary: false };
  }
  const upper = 1 - 1e-8;
  const qu = pearson(upper);
  if (qu >= degreesOfFreedom) {
    return {
      rho: upper,
      scale: Math.max(1, qu / degreesOfFreedom),
      pearson: qu,
      boundary: true,
    };
  }
  let low = 0;
  let high = upper;
  for (let iteration = 0; iteration < 100; iteration += 1) {
    const middle = (low + high) / 2;
    if (pearson(middle) > degreesOfFreedom) low = middle;
    else high = middle;
    if (high - low < 1e-10) break;
  }
  const rho = (low + high) / 2;
  return { rho, scale: 1, pearson: pearson(rho), boundary: false };
}

export function fitGuide(countInput, totalInput, design, options = {}) {
  const count = Array.from(countInput, Number);
  const total = Array.from(totalInput, Number);
  const maxIterations = options.maxIterations ?? 100;
  const tolerance = options.tolerance ?? 1e-8;
  const muBound = options.muBound ?? 1e-8;
  if (count.length !== total.length || count.length !== design.length) {
    throw new Error("Counts, totals, and design rows must have equal length.");
  }
  if (!count.length || !design[0]?.length) {
    throw new Error("A non-empty response and design matrix are required.");
  }
  const rank = design[0].length;
  const degreesOfFreedom = design.length - rank;
  if (degreesOfFreedom <= 0) {
    throw new Error("The model needs more samples than coefficients.");
  }
  for (let index = 0; index < count.length; index += 1) {
    if (!Number.isInteger(count[index]) || !Number.isInteger(total[index]) ||
        count[index] < 0 || total[index] <= 0 || count[index] > total[index]) {
      throw new Error("Counts must be integers satisfying 0 <= count <= total.");
    }
  }

  let beta = binomialInitial(count, total, design, tolerance);
  let converged = false;
  let iteration = 0;
  for (iteration = 1; iteration <= maxIterations; iteration += 1) {
    const eta = linearPredictor(design, beta);
    const mu = Float64Array.from(
      eta,
      (value) => clamp(logistic(value), muBound, 1 - muBound),
    );
    const rhoFit = estimateRho(count, total, mu, degreesOfFreedom);
    const response = new Float64Array(count.length);
    const weight = new Float64Array(count.length);
    for (let index = 0; index < count.length; index += 1) {
      response[index] = eta[index] +
        (count[index] / total[index] - mu[index]) /
        (mu[index] * (1 - mu[index]));
      weight[index] = total[index] * mu[index] * (1 - mu[index]) /
        (1 + (total[index] - 1) * rhoFit.rho);
    }
    const system = weightedSystem(design, weight, response);
    const next = solvePositiveDefinite(system.information, system.target);
    let change = 0;
    for (let index = 0; index < beta.length; index += 1) {
      change = Math.max(
        change,
        Math.abs(next[index] - beta[index]) / Math.max(1, Math.abs(beta[index])),
      );
    }
    beta = next;
    if (change < tolerance) {
      converged = true;
      break;
    }
  }

  const eta = linearPredictor(design, beta);
  const mu = Float64Array.from(
    eta,
    (value) => clamp(logistic(value), muBound, 1 - muBound),
  );
  const rhoFit = estimateRho(count, total, mu, degreesOfFreedom);
  const weight = Float64Array.from(
    mu,
    (value, index) => total[index] * value * (1 - value) /
      (1 + (total[index] - 1) * rhoFit.rho),
  );
  const system = weightedSystem(design, weight, eta);
  const covariance = inversePositiveDefinite(system.information);
  for (let row = 0; row < covariance.length; row += 1) {
    for (let column = 0; column < covariance.length; column += 1) {
      covariance[row][column] *= rhoFit.scale;
    }
  }
  const standardError = Float64Array.from(
    beta,
    (_, index) => Math.sqrt(covariance[index][index]),
  );
  const statistic = Float64Array.from(
    beta,
    (value, index) => value / standardError[index],
  );
  const pValue = Float64Array.from(
    statistic,
    (value) => studentTwoSidedP(value, degreesOfFreedom),
  );
  return {
    coefficient: Array.from(beta),
    standardError: Array.from(standardError),
    statistic: Array.from(statistic),
    pValue: Array.from(pValue),
    covariance: covariance.map((row) => Array.from(row)),
    fitted: Array.from(mu),
    rho: rhoFit.rho,
    scale: rhoFit.scale,
    pearson: rhoFit.pearson,
    dispersionBoundary: rhoFit.boundary,
    degreesOfFreedom,
    converged,
    iterations: Math.min(iteration, maxIterations),
  };
}

function inferVariableType(rows, variable) {
  const values = rows.map((row) => row[variable]);
  return values.every((value) => value !== "" && Number.isFinite(Number(value)))
    ? "numeric"
    : "categorical";
}

function expandVariable(rows, variable, requestedType) {
  const type = requestedType || inferVariableType(rows, variable);
  if (type === "numeric") {
    const values = rows.map((row) => Number(row[variable]));
    if (values.some((value) => !Number.isFinite(value))) {
      throw new Error(`${variable} contains a non-numeric or missing value.`);
    }
    return { columns: [variable], values: [values], reference: null, type };
  }
  const levels = [...new Set(rows.map((row) => String(row[variable])))]
    .sort((left, right) => left.localeCompare(right));
  if (levels.length < 2) {
    throw new Error(`${variable} must contain at least two levels.`);
  }
  const reference = levels[0];
  return {
    columns: levels.slice(1).map((level) => `${variable}${level}`),
    values: levels.slice(1).map(
      (level) => rows.map((row) => Number(String(row[variable]) === level)),
    ),
    reference,
    type,
  };
}

export function buildDesign(rows, configuration) {
  if (!rows.length) throw new Error("Sample metadata is empty.");
  const predictors = [
    configuration.predictor,
    ...(configuration.covariates || []),
  ].filter((value, index, array) => value && array.indexOf(value) === index);
  if (!predictors.length) throw new Error("Choose a primary predictor.");
  for (const variable of predictors) {
    if (!(variable in rows[0])) {
      throw new Error(`Metadata does not contain ${variable}.`);
    }
  }
  const expanded = new Map();
  for (const variable of predictors) {
    expanded.set(
      variable,
      expandVariable(rows, variable, configuration.types?.[variable]),
    );
  }
  const columnNames = ["(Intercept)"];
  const columnValues = [rows.map(() => 1)];
  for (const variable of predictors) {
    const item = expanded.get(variable);
    columnNames.push(...item.columns);
    columnValues.push(...item.values);
  }
  for (const interaction of configuration.interactions || []) {
    const [leftName, rightName] = interaction;
    if (!expanded.has(leftName) || !expanded.has(rightName) ||
        leftName === rightName) continue;
    const left = expanded.get(leftName);
    const right = expanded.get(rightName);
    for (let leftIndex = 0; leftIndex < left.columns.length; leftIndex += 1) {
      for (let rightIndex = 0; rightIndex < right.columns.length; rightIndex += 1) {
        columnNames.push(
          `${left.columns[leftIndex]}:${right.columns[rightIndex]}`,
        );
        columnValues.push(
          rows.map(
            (_, row) => left.values[leftIndex][row] *
              right.values[rightIndex][row],
          ),
        );
      }
    }
  }
  const matrix = rows.map(
    (_, row) => Float64Array.from(columnValues, (column) => column[row]),
  );
  const testWeight = new Float64Array(rows.length).fill(1);
  const testResponse = new Float64Array(rows.length).fill(1);
  const information = weightedSystem(matrix, testWeight, testResponse).information;
  cholesky(information);
  if (matrix.length <= columnNames.length) {
    throw new Error("The model needs more samples than coefficients.");
  }
  return {
    matrix,
    columns: columnNames,
    variables: Object.fromEntries(expanded),
    degreesOfFreedom: matrix.length - columnNames.length,
  };
}

export function runScreen(payload, progress = () => {}) {
  const {
    counts,
    guide,
    gene,
    control,
    totals,
    design,
    termIndex,
    minTotalCount = 10,
    fitOptions = {},
  } = payload;
  if (!Array.isArray(counts) || !counts.length) {
    throw new Error("The count matrix has no guide rows.");
  }
  if (!(termIndex >= 0 && termIndex < design[0].length)) {
    throw new Error("Choose a valid model coefficient.");
  }
  const results = [];
  const every = Math.max(1, Math.floor(counts.length / 100));
  for (let row = 0; row < counts.length; row += 1) {
    const totalCount = counts[row].reduce((sum, value) => sum + value, 0);
    const base = {
      guide: guide[row],
      gene: gene?.[row] ?? "",
      control: Boolean(control?.[row]),
      estimate: NaN,
      std_error: NaN,
      t_value: NaN,
      df: design.length - design[0].length,
      p_value: NaN,
      rho: NaN,
      mean_cpm: counts[row].reduce(
        (sum, value, column) => sum + value / totals[column] * 1e6,
        0,
      ) / totals.length,
      converged: false,
      dispersion_boundary: false,
    };
    if (totalCount >= minTotalCount) {
      try {
        const fit = fitGuide(counts[row], totals, design, fitOptions);
        base.estimate = fit.coefficient[termIndex];
        base.std_error = fit.standardError[termIndex];
        base.t_value = fit.statistic[termIndex];
        base.df = fit.degreesOfFreedom;
        base.p_value = fit.pValue[termIndex];
        base.rho = fit.rho;
        base.converged = fit.converged;
        base.dispersion_boundary = fit.dispersionBoundary;
      } catch {
        // One sparse or singular guide should not abort the complete screen.
      }
    }
    results.push(base);
    if (row % every === 0 || row === counts.length - 1) {
      progress((row + 1) / counts.length);
    }
  }
  const adjusted = bhAdjust(results.map((row) => row.p_value));
  results.forEach((row, index) => {
    row.fdr = adjusted[index];
  });
  return results;
}

export function calibrateControls(results, options = {}) {
  const alpha = options.alpha ?? 0.05;
  const minControls = options.minControls ?? 20;
  const minScale = options.minScale ?? 1;
  const valid = results.filter(
    (row) => row.control && Number.isFinite(row.t_value) &&
      Number.isFinite(row.df),
  );
  if (valid.length < minControls) {
    throw new Error(`At least ${minControls} finite negative controls are required.`);
  }
  const degrees = [...new Set(valid.map((row) => row.df))];
  if (degrees.length !== 1 || !(degrees[0] > 0)) {
    throw new Error("Negative controls must share one positive residual df.");
  }
  const empirical = quantileType8(
    valid.map((row) => Math.abs(row.t_value)),
    1 - alpha,
  );
  let low = 0;
  let high = 20;
  const target = alpha;
  for (let iteration = 0; iteration < 100; iteration += 1) {
    const middle = (low + high) / 2;
    if (studentTwoSidedP(middle, degrees[0]) > target) low = middle;
    else high = middle;
  }
  const reference = (low + high) / 2;
  const scale = Math.max(minScale, empirical / reference);
  const calibrated = results.map((row) => {
    const output = { ...row };
    output.raw_std_error = row.std_error;
    output.raw_t_value = row.t_value;
    output.raw_p_value = row.p_value;
    output.raw_fdr = row.fdr;
    if (Number.isFinite(row.std_error)) output.std_error = row.std_error * scale;
    if (Number.isFinite(row.t_value)) output.t_value = row.t_value / scale;
    if (Number.isFinite(output.t_value)) {
      output.p_value = studentTwoSidedP(output.t_value, row.df);
    }
    return output;
  });
  const adjusted = bhAdjust(calibrated.map((row) => row.p_value));
  calibrated.forEach((row, index) => {
    row.fdr = adjusted[index];
  });
  return { results: calibrated, scale, alpha, controls: valid.length };
}

export function geneConsistency(results, options = {}) {
  const minGuides = options.minGuides ?? 3;
  const alpha = options.alpha ?? 0.05;
  const minControlGenes = options.minControlGenes ?? 10;
  const minScale = options.minScale ?? 1;
  const grouped = new Map();
  for (const row of results) {
    const standardError = Number.isFinite(row.raw_std_error)
      ? row.raw_std_error
      : row.std_error;
    if (!row.gene || !Number.isFinite(row.estimate) ||
        !(standardError > 0)) continue;
    if (!grouped.has(row.gene)) grouped.set(row.gene, []);
    grouped.get(row.gene).push({ ...row, standardError });
  }
  const genes = [];
  for (const [gene, guides] of grouped) {
    if (guides.length < minGuides) continue;
    if (guides.some((guide) => guide.control) &&
        guides.some((guide) => !guide.control)) {
      throw new Error(`${gene} mixes control and non-control guides.`);
    }
    const weights = guides.map((guide) => 1 / (guide.standardError ** 2));
    const weightSum = weights.reduce((sum, value) => sum + value, 0);
    const estimate = guides.reduce(
      (sum, guide, index) => sum + weights[index] * guide.estimate,
      0,
    ) / weightSum;
    const stdError = Math.sqrt(1 / weightSum);
    const nonzero = guides.filter((guide) => guide.estimate !== 0);
    const agreement = estimate === 0 || !nonzero.length
      ? NaN
      : nonzero.filter((guide) => Math.sign(guide.estimate) === Math.sign(estimate))
        .length / nonzero.length;
    genes.push({
      gene,
      n_guides: guides.length,
      estimate,
      std_error: stdError,
      raw_statistic: estimate / stdError,
      guide_direction_agreement: agreement,
      converged_fraction: guides.filter((guide) => guide.converged).length /
        guides.length,
      control_gene: guides.every((guide) => guide.control),
    });
  }
  if (genes.length < 2) {
    throw new Error("At least two genes need enough finite guide scores.");
  }
  const raw = genes.map((gene) => gene.raw_statistic);
  const globalCenter = median(raw);
  let globalScale = median(raw.map((value) => Math.abs(value - globalCenter))) /
    normalQuantile(0.75);
  if (!(globalScale > 0) || !Number.isFinite(globalScale)) globalScale = 1;
  const controls = genes.filter((gene) => gene.control_gene)
    .map((gene) => gene.raw_statistic);
  const enoughControls = controls.length >= minControlGenes;
  const center = enoughControls ? median(controls) : globalCenter;
  const controlScale = enoughControls
    ? quantileType8(
      controls.map((value) => Math.abs(value - center)),
      1 - alpha,
    ) / normalQuantile(1 - alpha / 2)
    : 0;
  const scale = Math.max(minScale, globalScale, controlScale);
  for (const gene of genes) {
    gene.null_center = center;
    gene.null_scale = scale;
    gene.statistic = (gene.raw_statistic - center) / scale;
    gene.p_value = normalTwoSidedP(gene.statistic);
  }
  const adjusted = bhAdjust(genes.map((gene) => gene.p_value));
  genes.forEach((gene, index) => {
    gene.fdr = adjusted[index];
  });
  return {
    results: genes,
    center,
    scale,
    controlGenes: controls.length,
    usedControlNull: enoughControls,
  };
}
