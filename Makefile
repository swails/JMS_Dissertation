PAPER=JMS_Dissertation

install: paper

clean:
	/bin/rm -f *.aux *.log *~ $(PAPER).dvi *.sty *.def $(PAPER).pdf \
	        *.ilg $(PAPER).lof $(PAPER).lot $(PAPER).blg $(PAPER).out \
	        $(PAPER).toc $(PAPER).ps

veryclean: clean
	/bin/rm -f *.bbl

paper: clean
	@for f in DEFS/*; do ln -s $$f; done
	@for f in STYS/*; do ln -s $$f; done
	-pdflatex $(PAPER).tex
	-pdflatex $(PAPER).tex
	-pdflatex $(PAPER).tex
	-bibtex $(PAPER)
	-dvipdfmx $(PAPER).dvi
	/bin/rm *.sty *.def

