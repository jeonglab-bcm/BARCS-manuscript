# Raw benchmark data

Raw downloaded datasets are intentionally excluded from Git history. Place the
following files under this directory before rerunning the corresponding
benchmarks.

## GSE70038

Download the processed count matrix:

```sh
curl -L \
  'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE70038/suppl/GSE70038_rawReadCounts.txt.gz' \
  -o data/raw/GSE70038_rawReadCounts.txt.gz
```

## Chronos/Tzelepis

See [`chronos/README.md`](chronos/README.md) for the Tzelepis supplementary
files, Chronos Figshare file identifiers, and DepMap sources.

## Waterbear/GSE242880

See [`waterbear/README.md`](waterbear/README.md) for direct download URLs and
SHA-256 checksums.
