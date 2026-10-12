# BARCS methodology tutorial

A self-study guide to the BARCS 0.2.0 methodology, with worked examples, a
32-question multiple-choice quiz, and an answer key with explanations.

| File | Language | Build |
|---|---|---|
| `student_guide/barcs_student_guide.pdf` | English, grades 7-9 | `student_guide/build.sh`; both chapters in one book; start here |
| `student_guide/chapter1/barcs_student_guide_ch1.pdf` | English, grades 7-9 | Chapter 1, the ideas: counts, proportions, log-odds, noise, the design matrix (16-question quiz) |
| `student_guide/slides/barcs_slides_ch1.pdf`, `barcs_slides_ch2.pdf` | English, grades 7-9 | Beamer slides (16:9) that follow the two chapters step by step, with the same examples and numbers; no quiz |
| `student_guide/chapter2/barcs_student_guide_ch2.pdf` | English, grades 7-9 | Chapter 2, the calculations: fitting, noise, standard error, moderation, calibration, genes, FDR (12-question quiz) |
| `barcs_methods_tutorial.tex` / `.pdf` | English | `pdflatex` (or `latexmk -pdf`) |
| `barcs_methods_tutorial_ko.tex` / `.pdf` | Korean | `xelatex` (or `latexmk -xelatex`); needs `kotex` and the Nanum fonts |

The student guide is the introductory version. It follows ASD-STE100
Simplified Technical English (short sentences, active voice, every symbol
defined), carries one small data set through every step with all arithmetic
shown, and has a quiz at the end of each chapter whose answers show the calculation.
Each chapter builds as its own booklet; Chapter 2 opens with a one-page
recap of Chapter 1 so it can be handed out alone. The chapter text lives in
`chapter1/chapter1.tex` and `chapter2/chapter2.tex`; edit those, then run
`build.sh`.
Its numbers come from BARCS 0.2.1 (`bbreg()` on the four-library example).

The two methods tutorials use the same option order, so the answer keys are identical.
Formulas follow the BARCS package code (`R/bbreg.R`, `R/calibrate.R`,
`R/screen.R`, `R/gene.R`, `R/utils.R`).
