\name{NELS}
\alias{NELS88}
\title{National Education Longitudinal Study Data}
\description{
  Random sample of size 1000 from the US National Education Longitudinal
  Study (NELS) data containing the mathematics, science and reading scores,
  together with covariates, of 8th graders in 1988.
}
\usage{
data("NELS88")
}
\format{
  \code{\link{data.frame}} containing the identification
  number of the school to which the student belongs (\code{ID}),
  the standardized score of the student on a
  mathematics achievement test (\code{Math}; rescaled by an Item
  Response Theory (IRT) method where a higher score indicates greater
  proficiency in mathematics), the standardized
  score of the student on a science achievement test (\code{Science}),
  the standardized score of the student on a reading achievement test
  (\code{Reading}), a factor indicating whether the student is a member
  of an ethnic minority group (\code{Minority}), a numeric measure of
  the socio-economic status of the student and family (\code{SES}), a
  factor indicating whether the student is female (\code{Female}), a
  factor indicating whether the school is publicly funded
  (\code{Public}), the size of the student's school (\code{Size}),
  a factor indicating whether the school is located in an urban environment
  (\code{Urban}) and a factor indicating whether the school is located in a
  rural environment (\code{Rural}).
}
\source{
  Edward W. Frees, \sQuote{Student Achievement Data} in
  \url{https://sites.google.com/a/wisc.edu/jed-frees/multivariate-regression-using-copulas}.

  Originally, the National Center for Education Statistics page,
  \url{https://nces.ed.gov/surveys/nels88/}
}
%% \references{%% publications
%% }
\examples{
data("NELS88")
str(NELS88)
ftable(xtabs(~ Urban+Rural + Public, NELS88))#
## Add more sensible variable, ordered factor rural < agglo < urban
NELS88. <- within(NELS88, {
       UR <- factor(Urban:Rural, labels = c("agglo", "rural", "urban"))
      Urbanity <- ordered(UR, levels = c("rural", "agglo", "urban"))
      rm(UR) })
unique(NELS88.[, c("Urban","Rural", "Urbanity")]) # indeed, just 3 combination cases

xtabs(~ Minority+Urbanity, NELS88.) # (_not_ independent)
ftable(xtabs(~ Public+Urbanity+Female+Minority, NELS88.) -> tab.)
summary(tab.) # very very clearly not independent
}
\keyword{datasets}
