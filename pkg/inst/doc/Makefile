###--- Makefile mainly from Douglas Bates

## need to get 'texclean' command

DRIVERS = nacopula
#R_LIB = /usr/local/lib/R
R_LIB = /usr/local/app/R/current
R_PROG = R-devel
R_PROG = R
R_CALL = ${R_PROG} --vanilla --silent
#                  ^^^^^^^^^ also loses R_LIBS
R_CALL = ${R_PROG} --no-init-file --no-save --silent
R_CALL = ${R_PROG} --no-init-file --no-save --silent

TPUTS = ${TEXINPUTS}:${R_LIB}/share/texmf//

all: $(DRIVERS:=.pdf)
	${RM} *.out

%.tex: %.Rnw
	${R_PROG} CMD Sweave $<
##	echo "library(tools); Sweave('$<')" | ${R_CALL}
## I cannot yet set '--vanilla' when using 'R CMD Sweave' ...

%.R: %.Rnw
	${R_PROG} CMD Stangle $<

%.pdf : %.tex
	TEXINPUTS=${TPUTS} texi2dvi --pdf $<

%.dvi : %.tex
	TEXINPUTS=${TPUTS} texi2dvi $<


purge:
	$(RM) Rplots.pdf $(DRIVERS:=.out) *.aux *.log *.blg
clean: purge
	$(RM) $(DRIVERS:=.tex) $(DRIVERS:=.bbl)

realclean: clean
	$(RM) $(DRIVERS:=.pdf)
