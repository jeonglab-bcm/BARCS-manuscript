#!/usr/bin/env bash
# Build EPUB editions of the Korean BARCS methods tutorial.
#
#   teaching/epub/build_epub.sh
#
# Writes, next to the PDF in teaching/:
#   barcs_methods_tutorial_ko.epub         formulas as SVG images (any EPUB3 reader)
#   barcs_methods_tutorial_ko_mathml.epub  formulas as MathML (Apple Books, Kobo, Thorium)
# Needs pandoc >= 3, latex, dvisvgm and the Nanum fonts.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
teaching="$(dirname "$here")"
build="$(mktemp -d)"
trap 'rm -rf "$build"' EXIT

# A quick PDF pass gives the .aux with the printed equation numbers.
(cd "$build" && xelatex -no-pdf -interaction=batchmode "$teaching/barcs_methods_tutorial_ko.tex" > /dev/null 2>&1 || true)
python3 "$here/preprocess.py" "$teaching/barcs_methods_tutorial_ko.tex" "$build/tutorial.tex" \
  "$build/barcs_methods_tutorial_ko.aux"

fonts=/usr/share/fonts/truetype/nanum
common=(
  --from=latex --to=epub3 --toc --toc-depth=2 --split-level=1
  --metadata=title:"BARCS 이해하기"
  --metadata=subtitle:"방법론 자습 가이드: 예제, 객관식 퀴즈, 해설 (BARCS 버전 0.2.0)"
  --metadata=author:"정현환" --metadata=author:"이규원"
  --metadata=lang:ko-KR --metadata=date:2026-10-06
  --metadata=toc-title:"차례"
  --css="$here/epub.css"
  --epub-embed-font="$fonts/NanumMyeongjo.ttf"
  --epub-embed-font="$fonts/NanumMyeongjoBold.ttf"
)

(cd "$build" && MATH_DIR=math pandoc tutorial.tex "${common[@]}" \
   --lua-filter="$here/math-svg.lua" -o "$teaching/barcs_methods_tutorial_ko.epub")
(cd "$build" && pandoc tutorial.tex "${common[@]}" --mathml \
   -o "$teaching/barcs_methods_tutorial_ko_mathml.epub")
ls -l "$teaching"/barcs_methods_tutorial_ko*.epub
