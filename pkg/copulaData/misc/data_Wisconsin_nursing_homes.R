## Data preparation for the Wisconsin nursing homes utilization data of
## Sun, Frees, Rosenberg (2008, "Heavy-tailed longitudinal data modeling using copulas")


library(plyr)


### 1 Cleaning #################################################################

## Read the (raw) data from J. Sun
raw <- read.table(bzfile("data_Wisconsin_nursing_homes_raw.dat.bz2"),
                  sep = "", header = TRUE, na.string = ".")
stopifnot(length(unique(raw$POPID)) == 402)
dat <- subset(raw, !(POPID %in% c(577, 222, 234)))
dat <- subset(dat, CRYEAR >= 1995 & CRYEAR <= 2002)
dat <- subset(dat, F39 > 0)
stopifnot(length(unique(dat$POPID)) == 377) # now as in paper

## Renaming of some of the variables
dat <- within(dat, {
    ID        <- factor(POPID)
    Rate      <- Rate1
    LnNumBed  <- LNF681
    LnSqrFoot <- LNF39
    CRYear    <- factor(CRYEAR, ordered = TRUE) # => Year = CRYear - 1994
    TaxExempt <- factor(NP)
    SelfIns   <- factor(F23)
    MCert     <- factor(F26)
    Pro <- factor(Pro)
    Urban <- factor(Urban)
})

## Grab out those we work with
dat <- dat[, c("ID", "Rate", "LnNumBed", "LnSqrFoot", "CRYear",
               "Pro", "TaxExempt", "SelfIns", "MCert", "Urban")]
str(dat)


### 2 Save the cleaned data ####################################################

nursingHomes <- dat
str(nursingHomes)
save(nursingHomes, file = "nursingHomes.rda", compress = "xz")


### 3 Reproducing Table 1 and 2 of Sun et al. (2008) ###########################

## Table 1
tab1 <- ddply(dat, .(CRYear), summarize,
              Number  = length(Rate),
              Mean    = mean(Rate),
              Median  = median(Rate),
              Sd      = sd(Rate),
              Minimum = min(Rate),
              Maximum = max(Rate))
round(tab1, 2) # only the maxima differ a bit

## Table 2 entries
exp(median(dat$LnNumBed))
round(exp(median(dat$LnSqrFoot)), 2)
ddply(dat, .(SelfIns), summarize,
      Proportion = length(Rate) / nrow(dat),
      Median     = median(Rate))
ddply(dat, .(MCert), summarize,
      Proportion = length(Rate) / nrow(dat),
      Median     = median(Rate))
ddply(dat, .(Pro), summarize,
      Proportion = length(Rate) / nrow(dat),
      Median     = median(Rate))
ddply(dat, .(TaxExempt), summarize,
      Proportion = length(Rate) / nrow(dat),
      Median     = median(Rate))
ddply(dat, .(Urban), summarize,
      Proportion = length(Rate) / nrow(dat),
      Median     = median(Rate))


