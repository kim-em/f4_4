#!/usr/bin/env bash
# Builds the paper from the Verso source: TeX, BibTeX, PDF, checks, and the committed copies.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
job=f4-fsymbols
# Reproducible PDF: fixed dates and identifiers (1 November 2026, the month of the paper's date).
export SOURCE_DATE_EPOCH=${SOURCE_DATE_EPOCH:-1793491200}
export FORCE_SOURCE_DATE=1
lake build
rm -rf _out && mkdir -p _out/bib _out/pdf
# 1. TeX that reads the generated .bib, to let BibTeX format the cited references.
.lake/build/bin/paper --out _out/bib
(cd _out/bib && pdflatex -interaction=nonstopmode -halt-on-error $job.tex > pdflatex.log && bibtex $job > bibtex.log)
# 2. Self-contained TeX with the formatted bibliography included, compiled from scratch.
.lake/build/bin/paper --bbl _out/bib/$job.bbl --out _out/pdf
for _ in 1 2 3; do
  (cd _out/pdf && pdflatex -interaction=nonstopmode -halt-on-error $job.tex > pdflatex.log)
done
python3 check.py
cp _out/pdf/$job.tex ../$job.tex
cp _out/pdf/$job.pdf ../$job.pdf
