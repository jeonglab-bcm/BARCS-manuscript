#!/usr/bin/env python3
"""Apply MAGeCK's official piecewise CNV normalizer to external gene effects."""

from __future__ import annotations

import argparse
import csv

import numpy as np


_genfromtxt = np.genfromtxt


def _legacy_genfromtxt(*args, **kwargs):
    kwargs.setdefault("encoding", "bytes")
    return _genfromtxt(*args, **kwargs)


np.genfromtxt = _legacy_genfromtxt
np.float = float

from mageck.cnv_normalization import read_CNVdata, sgRNAscore_piecewisenorm


parser = argparse.ArgumentParser()
parser.add_argument("--effects", required=True)
parser.add_argument("--cnv", required=True)
parser.add_argument("--cell-line", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

genes: list[str] = []
effects: list[float] = []
with open(args.effects, newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        genes.append(row["gene"])
        effects.append(float(row["estimate"]))

cnv_array, _, cnv_gene = read_CNVdata(args.cnv, [args.cell_line])
corrected = sgRNAscore_piecewisenorm(
    effects, genes, cnv_array, cnv_gene, selectGenes=False
)

with open(args.output, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["gene", "cnv_corrected_estimate"])
    writer.writerows(zip(genes, corrected))
