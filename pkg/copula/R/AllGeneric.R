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

setGeneric("paramNames", function(x) standardGeneric("paramNames"))
setMethod("paramNames", "xcopula", function(x) paramNames(x@copula))

## get parameters
setGeneric("getParam", function(copula, freeOnly = TRUE, attr = FALSE, named = attr)
    standardGeneric("getParam"), signature = "copula")
## assign free parameters
setGeneric("freeParam<-", function(copula, value) standardGeneric("freeParam<-"))
## set or modify "fixedness" of parameters
setGeneric("fixedParam<-", function(copula, value) standardGeneric("fixedParam<-"))
## logical vector indicating which parameters are free
setGeneric("isFree", function(copula) standardGeneric("isFree"))

## number of (all/free) parameters
setGeneric("nParam", function(copula, freeOnly = FALSE) standardGeneric("nParam"),
           signature = "copula")


