export function formatNumber(value, digits = 3) {
  if (!Number.isFinite(value)) return "—";

  // Table renderers receive both (value, row). When this formatter is passed
  // directly, the row becomes the second argument. Treat anything other than
  // an integer precision as the default instead of forwarding it to Intl.
  const fractionDigits = Number.isInteger(digits)
    ? Math.max(0, Math.min(20, digits))
    : 3;
  const absolute = Math.abs(value);
  if (absolute !== 0 && (absolute < 0.001 || absolute >= 10000)) {
    return value.toExponential(2);
  }
  return value.toLocaleString(undefined, {
    maximumFractionDigits: fractionDigits,
  });
}

export function formatP(value) {
  if (!Number.isFinite(value)) return "—";
  if (value < 0.001) return value.toExponential(2);
  return value.toFixed(3);
}
