import {
  addReadToQuantification,
  buildLibraryIndex,
  createQuantification,
  determineLibrary,
  finalizeQuantification,
  stripFastqExtension,
} from "./fastq-core.js";

async function textStream(file) {
  let stream = file.stream();
  if (/\.gz$/i.test(file.name)) {
    if (typeof DecompressionStream === "undefined") {
      throw new Error("This browser cannot decompress gzip FASTQ files locally.");
    }
    stream = stream.pipeThrough(new DecompressionStream("gzip"));
  }
  return stream.pipeThrough(new TextDecoderStream());
}

async function visitFastq(file, visitor, maximumReads = Infinity) {
  const reader = (await textStream(file)).getReader();
  let buffer = "";
  let lines = [];
  let reads = 0;
  let done = false;
  while (!done && reads < maximumReads) {
    const chunk = await reader.read();
    done = chunk.done;
    buffer += chunk.value || "";
    const parts = buffer.split(/\r?\n/);
    buffer = parts.pop() || "";
    for (const part of parts) lines.push(part);
    while (lines.length >= 4 && reads < maximumReads) {
      const record = lines.splice(0, 4);
      if (!record[0].startsWith("@") || !record[2].startsWith("+") ||
          record[1].length !== record[3].length) {
        throw new Error(
          `${file.name} has an invalid FASTQ record near read ${reads + 1}.`,
        );
      }
      visitor(record[1], reads);
      reads += 1;
    }
  }
  if (reads >= maximumReads) {
    await reader.cancel();
  } else {
    if (buffer) lines.push(buffer);
    if (lines.some((line) => line !== "")) {
      throw new Error(`${file.name} ends with an incomplete FASTQ record.`);
    }
  }
  return reads;
}

async function sampleReads(file, maximumReads) {
  const reads = [];
  await visitFastq(file, (sequence) => reads.push(sequence), maximumReads);
  return reads;
}

async function quantifyFile(file, library, fileIndex, fileCount) {
  const index = buildLibraryIndex(library);
  const result = createQuantification(library);
  await visitFastq(file, (sequence, readIndex) => {
    addReadToQuantification(sequence, index, result);
    if (readIndex > 0 && readIndex % 50000 === 0) {
      self.postMessage({
        type: "progress",
        phase: "quantify",
        sample: stripFastqExtension(file.name),
        reads: readIndex,
        fraction: (fileIndex + 0.5) / fileCount,
      });
    }
  });
  return {
    name: stripFastqExtension(file.name),
    fileName: file.name,
    ...finalizeQuantification(result),
  };
}

self.onmessage = async ({ data }) => {
  if (data.type !== "quantify") return;
  try {
    const files = data.files;
    const libraries = data.libraries;
    if (!files.length || !libraries.length) {
      throw new Error("Choose FASTQ files and at least one guide library.");
    }
    self.postMessage({ type: "progress", phase: "detect", fraction: 0.01 });
    const reads = await sampleReads(files[0], data.maximumDetectionReads || 20000);
    const detection = determineLibrary(reads, libraries);
    const library = libraries.find((candidate) =>
      candidate.name === detection.selected
    );
    const samples = [];
    for (let index = 0; index < files.length; index += 1) {
      self.postMessage({
        type: "progress",
        phase: "quantify",
        sample: stripFastqExtension(files[index].name),
        reads: 0,
        fraction: index / files.length,
      });
      samples.push(await quantifyFile(files[index], library, index, files.length));
    }
    self.postMessage({
      type: "complete",
      result: { library, detection, samples },
    });
  } catch (error) {
    self.postMessage({
      type: "error",
      message: error instanceof Error ? error.message : String(error),
    });
  }
};
