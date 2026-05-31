PDF      := article.pdf
TEX      := article.tex
BIB      := references.bib
ENGINE   ?= tectonic

.PHONY: all clean view

all: $(PDF)

$(PDF): $(TEX) $(BIB) $(wildcard figures/*.pdf)
	$(ENGINE) $(TEX)

view: $(PDF)
	@command -v open >/dev/null && open $(PDF) || xdg-open $(PDF)

clean:
	rm -f article.aux article.bbl article.blg article.log article.out article.toc article.fls article.fdb_latexmk
	rm -f $(PDF)
