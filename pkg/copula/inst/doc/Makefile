###--- Makefile (adapted from pcalg & robustbase; origin from Douglas Bates)

## need to get 'texclean' command

DRIVERS = nacopula-pkg Frank-Rmpfr

#R_PROG = R-devel
R_PROG = R
#R_CALL = ${R_PROG} --vanilla --silent
#                  ^^^^^^^^^ also loses R_LIBS
R_CALL = ${R_PROG} --no-init-file --no-save --silent

## "Writing R Extensions" requires: 1st target must also clean up
## NB: This is *BAD* for submitting to JSS (which needs  $(DRIVERS).tex !)
all: $(DRIVERS:=.pdf)
	${MAKE} clean

paper-etc: $(DRIVERS:=.zip)
$(DRIVERS:=.zip): $(DRIVERS:=.Rnw) $(DRIVERS:=.pdf) $(DRIVERS:=.R) *.bib
	zip $@ $^
	ls -l $@ ; zipinfo $@

paper-all.zip: $(DRIVERS:=.Rnw) $(DRIVERS:=.tex) $(DRIVERS:=.pdf) $(DRIVERS:=.R) *.bib $(wildcard $(DRIVERS)-*.pdf)
	zip $@ $^
	ls -l $@ ; zipinfo $@

%.tex: %.Rnw
	${R_PROG} CMD Sweave $<

%.R: %.Rnw
	${R_PROG} CMD Stangle $<

%.pdf : %.tex *.bib
	${R_PROG} CMD texi2pdf $<
	qpdf --stream-data=compress $@ temp$@
	mv temp$@ $@

%.dvi : %.tex *.bib
	${R_PROG} CMD texi2dvi $<
## To make the paper (size in bytes) much smaller, use
##  make nacopula-pkg-ex2-splom.png ; pdflatex nacopula-pkg
%.png : %.pdf
	convert -quality 95 $< $@

purge:
	$(RM) Rplots.pdf $(DRIVERS:=.out) *.aux *.log *.blg
clean: purge
	$(RM) $(DRIVERS:=.tex) $(DRIVERS:=.bbl) $(DRIVERS:=-*.pdf)

realclean: clean
	$(RM) $(DRIVERS:=.pdf)
