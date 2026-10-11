#!/usr/bin/env bash
# Build the combined student guide and the two chapter booklets.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
build() {  # build <folder> <file>
  (cd "$here/$1" && for pass in 1 2; do
     pdflatex -interaction=nonstopmode -halt-on-error "$2.tex" > /dev/null || {
       grep -n '^!' -A5 "$2.log"; exit 1; }
   done && rm -f "$2".{aux,log,out,toc})
  echo "built $1/$2.pdf ($(pdfinfo "$here/$1/$2.pdf" | awk '/Pages/{print $2}') pages)"
}
build . barcs_student_guide
build chapter1 barcs_student_guide_ch1
build chapter2 barcs_student_guide_ch2
