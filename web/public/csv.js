function detectDelimiter(text) {
  const firstLine = text.split(/\r?\n/, 1)[0] || "";
  const commas = (firstLine.match(/,/g) || []).length;
  const tabs = (firstLine.match(/\t/g) || []).length;
  return tabs > commas ? "\t" : ",";
}

export function parseDelimited(text, delimiter = detectDelimiter(text)) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let index = 0; index < text.length; index += 1) {
    const character = text[index];
    if (quoted) {
      if (character === '"') {
        if (text[index + 1] === '"') {
          field += '"';
          index += 1;
        } else {
          quoted = false;
        }
      } else {
        field += character;
      }
      continue;
    }
    if (character === '"') {
      quoted = true;
    } else if (character === delimiter) {
      row.push(field);
      field = "";
    } else if (character === "\n") {
      row.push(field.replace(/\r$/, ""));
      if (row.some((value) => value !== "")) rows.push(row);
      row = [];
      field = "";
    } else {
      field += character;
    }
  }
  row.push(field.replace(/\r$/, ""));
  if (row.some((value) => value !== "")) rows.push(row);
  if (quoted) throw new Error("A quoted field is not closed.");
  if (!rows.length) throw new Error("The file is empty.");
  const width = rows[0].length;
  if (width < 2) throw new Error("The file must have at least two columns.");
  if (rows.some((value) => value.length !== width)) {
    throw new Error("Rows have different numbers of columns.");
  }
  return rows;
}

export function parseObjects(text) {
  const rows = parseDelimited(text);
  const header = rows[0].map((value) => value.trim());
  if (header.some((value) => !value)) throw new Error("Column names cannot be empty.");
  if (new Set(header).size !== header.length) {
    throw new Error("Column names must be unique.");
  }
  return {
    header,
    rows: rows.slice(1).map((values) => Object.fromEntries(
      header.map((name, index) => [name, values[index].trim()]),
    )),
  };
}

export function parseBarcsInputs(countText, metadataText) {
  const countsFile = parseObjects(countText);
  const metadataFile = parseObjects(metadataText);
  const guideColumn = ["guide", "sgRNA", "sgrna", "id"].find(
    (name) => countsFile.header.includes(name),
  );
  const geneColumn = ["gene", "Gene", "target"].find(
    (name) => countsFile.header.includes(name),
  );
  const controlColumn = ["control", "negative_control", "ntc"].find(
    (name) => countsFile.header.includes(name),
  );
  const sampleColumn = ["sample", "samples", "sample_name"].find(
    (name) => metadataFile.header.includes(name),
  );
  if (!guideColumn) {
    throw new Error("The count file needs a `guide`, `sgRNA`, or `id` column.");
  }
  if (!sampleColumn) {
    throw new Error(
      "The metadata file needs a `sample`, `samples`, or `sample_name` column.",
    );
  }
  if (!countsFile.rows.length || !metadataFile.rows.length) {
    throw new Error("Both files need at least one data row.");
  }
  const samples = metadataFile.rows.map((row) => row[sampleColumn]);
  if (samples.some((sample) => !sample) ||
      new Set(samples).size !== samples.length) {
    throw new Error("Metadata sample names must be non-empty and unique.");
  }
  const missingSamples = samples.filter(
    (sample) => !countsFile.header.includes(sample),
  );
  if (missingSamples.length) {
    throw new Error(
      `Count columns are missing metadata samples: ${missingSamples.join(", ")}.`,
    );
  }
  const guides = countsFile.rows.map((row) => row[guideColumn]);
  if (guides.some((guide) => !guide) ||
      new Set(guides).size !== guides.length) {
    throw new Error("Guide names must be non-empty and unique.");
  }
  const counts = countsFile.rows.map((row, rowIndex) => samples.map((sample) => {
    const value = Number(row[sample]);
    if (!Number.isInteger(value) || value < 0) {
      throw new Error(
        `Count for ${guides[rowIndex]} in ${sample} is not a non-negative integer.`,
      );
    }
    return value;
  }));
  const totalColumn = ["total", "library_total", "library_size"].find(
    (name) => metadataFile.header.includes(name),
  );
  const totals = totalColumn
    ? metadataFile.rows.map((row) => Number(row[totalColumn]))
    : samples.map(
      (_, column) => counts.reduce((sum, row) => sum + row[column], 0),
    );
  if (totals.some((total) => !Number.isInteger(total))) {
    throw new Error(
      `Metadata ${totalColumn} values must be integer full-library totals.`,
    );
  }
  if (totals.some((total) => total <= 0)) {
    throw new Error("Every sample must have a positive full-library total.");
  }
  counts.forEach((row, rowIndex) => row.forEach((value, column) => {
    if (value > totals[column]) {
      throw new Error(
        `${guides[rowIndex]} exceeds the full-library total for ${samples[column]}.`,
      );
    }
  }));
  const truthy = new Set(["1", "true", "yes", "y", "control", "ntc"]);
  return {
    counts,
    totals,
    guide: guides,
    gene: countsFile.rows.map((row) => geneColumn ? row[geneColumn] : ""),
    control: countsFile.rows.map(
      (row) => truthy.has(
        String(controlColumn ? row[controlColumn] : "").toLowerCase(),
      ),
    ),
    samples,
    metadata: metadataFile.rows,
    metadataColumns: metadataFile.header.filter(
      (name) => name !== sampleColumn && name !== totalColumn,
    ),
    countColumns: countsFile.header,
    totalsSource: totalColumn || "column sums",
  };
}

function escapeField(value) {
  const text = value == null || Number.isNaN(value) ? "" : String(value);
  return /[",\n\r]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

export function toCsv(rows, preferredColumns = []) {
  if (!rows.length) return "";
  const discoveredSet = new Set();
  for (const row of rows) {
    for (const key of Object.keys(row)) discoveredSet.add(key);
  }
  const discovered = [...discoveredSet];
  const columns = [
    ...preferredColumns.filter((column) => discovered.includes(column)),
    ...discovered.filter((column) => !preferredColumns.includes(column)),
  ];
  const lines = [columns.map(escapeField).join(",")];
  for (const row of rows) {
    lines.push(
      columns.map((column) => escapeField(row[column])).join(","),
    );
  }
  return lines.join("\n");
}
