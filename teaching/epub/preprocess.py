"""Rewrite the tutorial's LaTeX into constructs pandoc understands.

    python3 preprocess.py <in.tex> <out.tex> [<in.aux>]

The PDF uses tcolorbox environments, an enumitem list with (A)-(D) labels and
kotex font commands. Pandoc ignores or drops these, so they are mapped to
block quotes with a bold title and to enumerate lists with explicit labels.
EPUB has no equation anchors, so with the PDF build's .aux the equation
labels become \tag{N} and \eqref{...} becomes the printed number "(N)".
"""
import re
import sys

source, target = sys.argv[1], sys.argv[2]
aux = sys.argv[3] if len(sys.argv) > 3 else None
text = open(source, encoding="utf-8").read()


def braced(s, start):
    """Return (content, end) of the {...} group opening at s[start]."""
    depth = 0
    for i in range(start, len(s)):
        depth += {"{": 1, "}": -1}.get(s[i], 0)
        if depth == 0:
            return s[start + 1:i], i + 1
    raise ValueError("unbalanced braces")


# \tfrac becomes an mfrac attribute the EPUB MathML schema rejects.
while "\\tfrac{" in text:
    i = text.index("\\tfrac{")
    num, j = braced(text, i + 6)
    den, k = braced(text, j)
    text = text[:i] + "{\\textstyle\\frac{" + num + "}{" + den + "}}" + text[k:]

# Equation numbers from the PDF build.
if aux:
    numbers = dict(re.findall(r"\\newlabel\{(eq:[^}]+)\}\{\{([^}]+)\}", open(aux, encoding="utf-8").read()))
    for label, number in numbers.items():
        text = text.replace(f"\\label{{{label}}}", f"\\tag{{{number}}}")
        text = text.replace(f"\\eqref{{{label}}}", f"({number})")

# Formulas in headings would put images into the table of contents.
text = re.sub(r"(\\(?:sub)*section\{[^\n]*)",
              lambda m: re.sub(r"\$([A-Za-z])\$", r"\\emph{\1}", m.group(1)), text)

# Font and box definitions only matter to the PDF.
text = re.sub(r"\\set(main|sans|mono)hangulfont\{[^}]*\}\n", "", text)
text = re.sub(r"\\newtcolorbox\{.*?\n(?:.*\n)*?.*?\}\n(?=\\|\n)", "", text)

titles = {"keyidea": "핵심 개념", "incode": "패키지에서는"}
for env, title in titles.items():
    text = text.replace(f"\\begin{{{env}}}", f"\\begin{{quote}}\\textbf{{{title}}}\\par ")
    text = text.replace(f"\\end{{{env}}}", "\\end{quote}")
# example takes an optional suffix: \begin{example}[ 2: ...]
text = re.sub(r"\\begin\{example\}(?:\[([^\]]*)\])?",
              lambda m: "\\begin{quote}\\textbf{예제" + (m.group(1) or "") + "}\\par ", text)
text = text.replace("\\end{example}", "\\end{quote}")

# Multiple-choice options: (A), (B), ...
text = text.replace("\\begin{choices}", "\\begin{enumerate}[(A)]")
text = text.replace("\\end{choices}", "\\end{enumerate}")
# key=value enumitem options confuse pandoc; keep the list, drop the options.
text = re.sub(r"\\begin\{enumerate\}\[label=[^\n]*\]", r"\\begin{enumerate}", text)

# texmath (MathML) needs braced arguments: \text c -> \text{c}.
text = re.sub(r"\\text ([A-Za-z0-9])", r"\\text{\1}", text)

# Title block: pandoc gets metadata from the build script instead.
text = re.sub(r"\\title\{.*?\}\n\\author\{.*?\}\n\\date\{.*?\}\n", "", text, flags=re.S)
text = text.replace("\\maketitle", "").replace("\\tableofcontents", "")
text = text.replace("\\newpage", "")
open(target, "w", encoding="utf-8").write(text)
