#!/usr/bin/env bash
#
# Install the pinned Julia toolchain into `.tools/julia`.
#
# Julia is the one dependency pixi cannot supply. conda-forge builds it for
# linux-64 only, and the CRISPulator simulations that produce every FACS
# benchmark input have to run on this machine, which is aarch64. So instead of
# a conda package we fetch the official binary for the host architecture and
# verify it against a checksum published by the Julia project.
#
# The version is pinned rather than floated because CRISPulator draws its
# screens from Julia's global RNG, and the exact stream is only guaranteed
# within a Julia release series. A different Julia gives different screens,
# which would silently invalidate every committed metric.
#
# Idempotent: re-running with the pinned version already installed does
# nothing. To upgrade, change JULIA_VERSION and the two checksums together --
# and expect the simulated inputs, and therefore the committed results, to
# change.

set -euo pipefail

JULIA_VERSION="1.12.6"

# From https://julialang-s3.julialang.org/bin/checksums/julia-${JULIA_VERSION}.sha256
SHA256_AARCH64="029b93b857bd0ffd627f9a8580d3bbaa63daf008d7b7aed02fbceb8fd57c4899"
SHA256_X86_64="bbabf3bef19421a9dbd24a767d807606ab85e444323b5a1c73ffe293fa3d079a"

# `pixi run` sets PIXI_PROJECT_ROOT. The fallback lets the script also be run
# directly from a shell that is not inside a pixi task.
ROOT="${PIXI_PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PREFIX="${ROOT}/.tools/julia"

case "$(uname -m)" in
  aarch64|arm64) JULIA_ARCH="aarch64"; JULIA_DIR="aarch64";  SHA256="${SHA256_AARCH64}" ;;
  x86_64)        JULIA_ARCH="x86_64";  JULIA_DIR="x86_64";   SHA256="${SHA256_X86_64}"  ;;
  *) echo "setup_julia.sh: unsupported architecture $(uname -m)" >&2; exit 1 ;;
esac

if [[ -x "${PREFIX}/bin/julia" ]] &&
   "${PREFIX}/bin/julia" --version 2>/dev/null | grep -q "${JULIA_VERSION}"; then
  echo "Julia ${JULIA_VERSION} already installed at ${PREFIX}"
  exit 0
fi

TARBALL="julia-${JULIA_VERSION}-linux-${JULIA_ARCH}.tar.gz"
URL="https://julialang-s3.julialang.org/bin/linux/${JULIA_DIR}/${JULIA_VERSION%.*}/${TARBALL}"

WORKDIR="$(mktemp -d)"
trap 'rm -rf "${WORKDIR}"' EXIT

echo "Downloading ${URL}"
curl -fsSL --retry 3 -o "${WORKDIR}/${TARBALL}" "${URL}"

echo "${SHA256}  ${WORKDIR}/${TARBALL}" | sha256sum --check --status || {
  echo "setup_julia.sh: checksum mismatch for ${TARBALL}; refusing to install" >&2
  exit 1
}

# Unpack to a temporary name and swap it in, so an interrupted extraction
# cannot leave a half-populated `.tools/julia` that the check above would then
# treat as present.
rm -rf "${PREFIX}" "${PREFIX}.incoming"
mkdir -p "${PREFIX}.incoming"
tar -xzf "${WORKDIR}/${TARBALL}" -C "${PREFIX}.incoming" --strip-components=1
mv "${PREFIX}.incoming" "${PREFIX}"

echo "Installed $("${PREFIX}/bin/julia" --version) at ${PREFIX}"
