### ../DESCRIPTION -- Collate:  want setClass() and setGeneric() first
### ~~~~~~~~~~~~~~    --------

## as we have ./Classes.R  ./AllClass.R  before this,
## currently, many setGeneric()s are there

##' A  (generic!) Function with methods to replace the 'fullname' slot :
setGeneric("describeCop", function(x, kind = c("short", "very short", "long"))
    standardGeneric("describeCop"))

## This could go back to ./dC-dc.R  if all "dCdu" methods also go there
##                       ~~~~~~~~~
setGeneric("dCdu", function(copula, u, ...) standardGeneric("dCdu"))


setGeneric("fitCopula", function(copula, data, ...) standardGeneric("fitCopula"))

setGeneric("gofCopula", function(copula, x, ...) standardGeneric("gofCopula"))

