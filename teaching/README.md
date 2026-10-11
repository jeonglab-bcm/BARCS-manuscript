# BARCS methodology tutorial

A self-study guide to the BARCS 0.2.0 methodology, with worked examples, a
32-question multiple-choice quiz, and an answer key with explanations.

| File | Language | Build |
|---|---|---|
| `barcs_methods_tutorial.tex` / `.pdf` | English | `pdflatex` (or `latexmk -pdf`) |
| `barcs_methods_tutorial_ko.tex` / `.pdf` | Korean | `xelatex` (or `latexmk -xelatex`); needs `kotex` and the Nanum fonts |
| `barcs_methods_tutorial_ko.epub` | Korean, e-reader | `epub/build_epub.sh`; formulas as SVG images, works on any EPUB3 reader |
| `barcs_methods_tutorial_ko_mathml.epub` | Korean, e-reader | same script; formulas as MathML (sharper on Apple Books, Kobo, Thorium; no equation numbers) |

Both versions use the same option order, so the answer keys are identical.
Formulas follow the BARCS package code (`R/bbreg.R`, `R/calibrate.R`,
`R/screen.R`, `R/gene.R`, `R/utils.R`).

The EPUB editions are generated from the Korean `.tex` with pandoc
(`epub/preprocess.py` maps the PDF-only boxes and lists to plain LaTeX,
`epub/math-svg.lua` renders formulas with xelatex and dvisvgm). Both pass
epubcheck 5.1 and embed NanumMyeongjo, so Korean text displays on readers
without a Korean font. Rebuild them after editing the `.tex`.
