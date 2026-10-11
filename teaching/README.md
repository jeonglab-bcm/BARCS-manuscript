# BARCS methodology tutorial

A self-study guide to the BARCS 0.2.0 methodology, with worked examples, a
32-question multiple-choice quiz, and an answer key with explanations.

| File | Language | Build |
|---|---|---|
| `barcs_student_guide.tex` / `.pdf` | English, grades 7-9 | `pdflatex` (run twice); start here |
| `barcs_methods_tutorial.tex` / `.pdf` | English | `pdflatex` (or `latexmk -pdf`) |
| `barcs_methods_tutorial_ko.tex` / `.pdf` | Korean | `xelatex` (or `latexmk -xelatex`); needs `kotex` and the Nanum fonts |

The student guide is the introductory version. It follows ASD-STE100
Simplified Technical English (short sentences, active voice, every symbol
defined), carries one small data set through every step with all arithmetic
shown, and ends with a 20-question quiz whose answers show the calculation.
Its numbers come from BARCS 0.2.1 (`bbreg()` on the four-library example).

The two methods tutorials use the same option order, so the answer keys are identical.
Formulas follow the BARCS package code (`R/bbreg.R`, `R/calibrate.R`,
`R/screen.R`, `R/gene.R`, `R/utils.R`).
