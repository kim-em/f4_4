all: build/f4-fsymbols.pdf

build/f4-fsymbols.pdf: f4-fsymbols.tex f4-fsymbols.bib
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build f4-fsymbols.tex

clean:
	latexmk -C -outdir=build f4-fsymbols.tex

.PHONY: all clean
