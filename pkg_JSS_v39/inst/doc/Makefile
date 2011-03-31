###--- Makefile mainly from Douglas Bates

## need to get 'texclean' command

DRIVERS = nacopula-pkg
#R_LIB = /usr/local/lib/R
R_LIB = /usr/local/app/R/current
#R_PROG = R-devel
R_PROG = R
#R_CALL = ${R_PROG} --vanilla --silent
#                  ^^^^^^^^^ also loses R_LIBS
R_CALL = ${R_PROG} --no-init-file --no-save --silent

TPUTS = ${TEXINPUTS}:${R_LIB}/share/texmf//

## "Writing R Extensions" requires: 1st target must also clean up
## NB: This is *BAD* for submitting to JSS (which needs  $(DRIVERS).tex !)
all: $(DRIVERS:=.pdf)
	${MAKE} clean

%.tex: %.Rnw
	${R_PROG} CMD Sweave $<

%.R: %.Rnw
	${R_PROG} CMD Stangle $<

%.pdf : %.tex *.bib
	TEXINPUTS=${TPUTS} texi2pdf $<

%.dvi : %.tex *.bib
	TEXINPUTS=${TPUTS} texi2dvi $<
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
