PAPER=JMS_Dissertation

install: paper

dissertation: paper

clean:
	/bin/rm -f *.aux *.log *~ $(PAPER).dvi *.sty *.def \
	        *.ilg $(PAPER).lof $(PAPER).lot $(PAPER).blg $(PAPER).out \
	        $(PAPER).toc $(PAPER).ps

veryclean: clean
	/bin/rm -f *.bbl $(PAPER).pdf $(PAPER).bib

paper: clean bib
	@for f in DEFS/*; do ln -s $$f; done
	@for f in STYS/*; do ln -s $$f; done
	-yes s | pdflatex $(PAPER).tex
	-yes s | pdflatex $(PAPER).tex
	-bibtex $(PAPER)
	-yes s | pdflatex $(PAPER).tex
	-dvipdfmx $(PAPER).dvi

bib::
	-/bin/rm -f $(PAPER).bib
	-cp $(GIT_HOME)/Bibliography/BibTeX/YorkLib.bib $(PAPER).bib
	-cat Extra.bib >> $(PAPER).bib
