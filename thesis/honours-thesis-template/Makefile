# Makefile for theses using acadia-masters-thesis.sty or acadia-hon-thesis.sty

COMPILE = pdflatex

BIB = bibtex

FILES = *.tex acadia*thesis.sty *.bib

AUX_FILES = *.aux *.toc *.bbl *.blg *.loa *.lof *.log *.lol *.lot *.out

PAPER = thesis

$(PAPER).pdf: $(FILES)
	$(COMPILE) $(PAPER)
	$(BIB) $(PAPER) || true
	$(COMPILE) $(PAPER)
	$(COMPILE) $(PAPER)



$(INSTRUCTION_PAPER).pdf: $(FILES)
	$(COMPILE) $(INSTRUCTION_PAPER)
	$(BIB) $(INSTRUCTION_PAPER) || true
	$(COMPILE) $(INSTRUCTION_PAPER)
	$(COMPILE) $(INSTRUCTION_PAPER)

all: $(PAPER).pdf

clean:
	rm -f $(AUX_FILES) $(PAPER).pdf
