import { parseObjects, toCsv } from "./csv.js";

const DNA = /^[ACGT]+$/i;
const GUIDE_ALIASES = ["guide", "sgRNA", "sgrna", "id"];
const SEQUENCE_ALIASES = ["sequence", "seq", "spacer"];
const GENE_ALIASES = ["gene", "Gene", "target"];
const CONTROL_ALIASES = ["control", "negative_control", "ntc"];
const BASE_CODE = Object.freeze({ A: 0, C: 1, G: 3, T: 2 });
const COMPLEMENT_CODE = Object.freeze([2, 3, 0, 1]);
const MAX_SAFE_KMER_LENGTH = 26;

function firstColumn(header, aliases) {
  return aliases.find((name) => header.includes(name));
}

function inferGene(guide) {
  const separator = guide.lastIndexOf("_");
  return separator > 0 ? guide.slice(0, separator) : guide;
}

function inferControl(guide, gene) {
  return /^(neg(?:ative)?(?:ctrl|control)?|ntc|non[-_ ]?target)/i.test(
    `${gene} ${guide}`,
  );
}

function logical(value, fallback = false) {
  if (value == null || value === "") return fallback;
  return new Set(["1", "true", "yes", "y", "control", "ntc"]).has(
    String(value).trim().toLowerCase(),
  );
}

function validateGuideRows(rows, name) {
  if (!rows.length) throw new Error(`${name} contains no guides.`);
  const ids = new Set();
  for (const row of rows) {
    row.guide = String(row.guide || "").trim();
    row.sequence = String(row.sequence || "").replace(/\s+/g, "").toUpperCase();
    row.gene = String(row.gene || inferGene(row.guide)).trim();
    row.control = Boolean(row.control);
    if (!row.guide || ids.has(row.guide)) {
      throw new Error(`${name} has an empty or repeated guide identifier.`);
    }
    if (!DNA.test(row.sequence) || row.sequence.length < 15 ||
        row.sequence.length > 40) {
      throw new Error(
        `${row.guide} must contain a 15–40 nt A/C/G/T guide sequence.`,
      );
    }
    ids.add(row.guide);
  }
  const frequencies = new Map();
  rows.forEach((row) => {
    frequencies.set(row.sequence, (frequencies.get(row.sequence) || 0) + 1);
  });
  const duplicateSequences = [...frequencies.values()]
    .filter((frequency) => frequency > 1)
    .reduce((sum, frequency) => sum + frequency, 0);
  const guides = rows.filter((row) => frequencies.get(row.sequence) === 1);
  if (!guides.length) {
    throw new Error(`${name} has no uniquely identifiable guide sequences.`);
  }
  const lengths = [...new Set(guides.map((row) => row.sequence.length))]
    .sort((left, right) => left - right);
  if (lengths.length !== 1) {
    throw new Error(
      `${name} contains mixed guide lengths (${lengths.join(", ")} nt). ` +
      "CB2-compatible k-mer counting requires one guide length per library.",
    );
  }
  return {
    name,
    guides,
    totalGuides: rows.length,
    duplicateSequences,
    excludedGuides: rows.length - guides.length,
    lengths,
  };
}

function parseFasta(text, name) {
  const rows = [];
  let header = null;
  let sequence = "";
  const finish = () => {
    if (!header) return;
    const fields = header.trim().split(/\s+/);
    const guide = fields.shift();
    const annotations = Object.fromEntries(fields.flatMap((field) => {
      const separator = field.indexOf("=");
      return separator > 0
        ? [[field.slice(0, separator).toLowerCase(), field.slice(separator + 1)]]
        : [];
    }));
    const gene = annotations.gene || annotations.target || inferGene(guide);
    rows.push({
      guide,
      sequence,
      gene,
      control: logical(
        annotations.control ?? annotations.ntc,
        inferControl(guide, gene),
      ),
    });
  };
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line.startsWith(">")) {
      finish();
      header = line.slice(1);
      sequence = "";
    } else {
      if (!header) throw new Error(`${name} is not valid FASTA.`);
      sequence += line;
    }
  }
  finish();
  return rows;
}

function parseGuideTable(text, name) {
  const table = parseObjects(text);
  const guideColumn = firstColumn(table.header, GUIDE_ALIASES);
  const sequenceColumn = firstColumn(table.header, SEQUENCE_ALIASES);
  const geneColumn = firstColumn(table.header, GENE_ALIASES);
  const controlColumn = firstColumn(table.header, CONTROL_ALIASES);
  if (!guideColumn || !sequenceColumn) {
    throw new Error(
      `${name} needs guide/id and sequence/spacer columns.`,
    );
  }
  return table.rows.map((row) => {
    const guide = row[guideColumn];
    const gene = geneColumn ? row[geneColumn] : inferGene(guide);
    return {
      guide,
      sequence: row[sequenceColumn],
      gene,
      control: controlColumn
        ? logical(row[controlColumn])
        : inferControl(guide, gene),
    };
  });
}

export function parseGuideLibrary(text, name = "guide library") {
  const trimmed = text.trim();
  const rows = trimmed.startsWith(">")
    ? parseFasta(trimmed, name)
    : parseGuideTable(trimmed, name);
  return validateGuideRows(rows, name);
}

export function reverseComplement(sequence) {
  const complement = { A: "T", C: "G", G: "C", T: "A", N: "N" };
  let output = "";
  for (let index = sequence.length - 1; index >= 0; index -= 1) {
    output += complement[sequence[index].toUpperCase()] || "N";
  }
  return output;
}

function createKmerCodec(length) {
  const big = length > MAX_SAFE_KMER_LENGTH;
  return {
    big,
    zero: big ? 0n : 0,
    radix: big ? 4n : 4,
    prefixModulus: big
      ? 4n ** BigInt(length - 1)
      : 4 ** (length - 1),
  };
}

function appendBase(code, base, observedLength, codec) {
  const value = codec.big ? BigInt(base) : base;
  if (observedLength < codec.length) {
    return code * codec.radix + value;
  }
  return (code % codec.prefixModulus) * codec.radix + value;
}

function encodeKmer(sequence, codec) {
  let code = codec.zero;
  for (const nucleotide of sequence) {
    code = code * codec.radix +
      (codec.big ? BigInt(BASE_CODE[nucleotide]) : BASE_CODE[nucleotide]);
  }
  return code;
}

export function buildLibraryIndex(library) {
  const length = library.lengths[0];
  const codec = { ...createKmerCodec(length), length };
  const guideByCode = new Map();
  library.guides.forEach((guide, guideIndex) => {
    guideByCode.set(encodeKmer(guide.sequence, codec), guideIndex);
  });
  return { library, length, codec, guideByCode };
}

function firstKmerHit(sequence, index, reverse = false) {
  const { codec, guideByCode, length } = index;
  let code = codec.zero;
  let observedLength = 0;
  for (let offset = 0; offset < sequence.length; offset += 1) {
    const sequenceIndex = reverse ? sequence.length - offset - 1 : offset;
    const base = BASE_CODE[sequence[sequenceIndex]];
    if (base === undefined) {
      code = codec.zero;
      observedLength = 0;
      continue;
    }
    const encodedBase = reverse ? COMPLEMENT_CODE[base] : base;
    code = appendBase(code, encodedBase, observedLength, codec);
    observedLength = Math.min(observedLength + 1, length);
    if (observedLength === length) {
      const guideIndex = guideByCode.get(code);
      if (guideIndex !== undefined) {
        return {
          guideIndex,
          orientation: reverse ? "-" : "+",
          position: reverse
            ? sequence.length - offset - 1
            : offset - length + 1,
        };
      }
    }
  }
  return null;
}

export function alignRead(read, index) {
  const sequence = read.trim().toUpperCase();
  const forward = firstKmerHit(sequence, index);
  if (forward) return { status: "mapped", ...forward };
  const reverse = firstKmerHit(sequence, index, true);
  return reverse
    ? { status: "mapped", ...reverse }
    : { status: "unmapped" };
}

function gini(values) {
  const sorted = [...values].sort((left, right) => left - right);
  const sum = sorted.reduce((total, value) => total + value, 0);
  if (!sum) return NaN;
  const weighted = sorted.reduce(
    (total, value, index) => total + (index + 1) * value,
    0,
  );
  return (2 * weighted) / (sorted.length * sum) -
    (sorted.length + 1) / sorted.length;
}

export function createQuantification(library) {
  return {
    counts: new Array(library.guides.length).fill(0),
    totalReads: 0,
    mappedReads: 0,
    ambiguousReads: 0,
    unmappedReads: 0,
    forwardReads: 0,
    reverseReads: 0,
    positionCounts: new Map(),
  };
}

export function addReadToQuantification(sequence, index, result) {
  result.totalReads += 1;
  const hit = alignRead(sequence, index);
  if (hit.status === "mapped") {
    result.counts[hit.guideIndex] += 1;
    result.mappedReads += 1;
    if (hit.orientation === "+") result.forwardReads += 1;
    else result.reverseReads += 1;
    result.positionCounts.set(
      hit.position,
      (result.positionCounts.get(hit.position) || 0) + 1,
    );
  } else {
    result.unmappedReads += 1;
  }
}

export function finalizeQuantification(result) {
  const positions = [...result.positionCounts.entries()]
    .sort((left, right) => right[1] - left[1]);
  const zeroGuides = result.counts.filter((count) => count === 0).length;
  return {
    ...result,
    positionCounts: positions,
    dominantPosition: positions[0]?.[0] ?? null,
    mappingRate: result.totalReads
      ? result.mappedReads / result.totalReads
      : 0,
    zeroGuides,
    zeroGuideFraction: zeroGuides / result.counts.length,
    gini: gini(result.counts),
  };
}

export function quantifyFastqText(text, library) {
  const lines = text.split(/\r?\n/);
  while (lines.length && lines.at(-1) === "") lines.pop();
  if (lines.length % 4 !== 0) {
    throw new Error("FASTQ must contain complete four-line records.");
  }
  const index = buildLibraryIndex(library);
  const result = createQuantification(library);
  for (let line = 0; line < lines.length; line += 4) {
    if (!lines[line].startsWith("@") || !lines[line + 2].startsWith("+")) {
      throw new Error(`Invalid FASTQ record beginning at line ${line + 1}.`);
    }
    if (lines[line + 1].length !== lines[line + 3].length) {
      throw new Error(`Sequence and quality lengths differ at line ${line + 1}.`);
    }
    addReadToQuantification(lines[line + 1], index, result);
  }
  return finalizeQuantification(result);
}

export function determineLibrary(reads, libraries) {
  if (!libraries.length) throw new Error("At least one guide library is required.");
  if (!reads.length) throw new Error("The FASTQ sample contains no reads.");
  const scores = libraries.map((library) => {
    const index = buildLibraryIndex(library);
    let mapped = 0;
    for (const read of reads) {
      const hit = alignRead(read, index);
      if (hit.status === "mapped") mapped += 1;
    }
    return {
      name: library.name,
      mapped,
      sampledReads: reads.length,
      mappingRate: reads.length ? mapped / reads.length : 0,
    };
  }).sort((left, right) =>
    right.mapped - left.mapped ||
    left.name.localeCompare(right.name)
  );
  if (scores[0].mapped === 0) {
    throw new Error(
      "No exact guide match was found in the sampled reads for any candidate library.",
    );
  }
  return {
    selected: scores[0].name,
    scores,
    confidence: scores.length < 2 || scores[1].mapped === 0
      ? Infinity
      : scores[0].mapped / scores[1].mapped,
  };
}

export function stripFastqExtension(filename) {
  return filename.replace(/\.(?:fastq|fq)(?:\.gz)?$/i, "");
}

export function quantificationToCountsText(library, samples) {
  return toCsv(
    library.guides.map((guide, guideIndex) => ({
      guide: guide.guide,
      gene: guide.gene,
      control: guide.control,
      ...Object.fromEntries(
        samples.map((sample) => [sample.name, sample.counts[guideIndex]]),
      ),
    })),
    ["guide", "gene", "control", ...samples.map((sample) => sample.name)],
  );
}
