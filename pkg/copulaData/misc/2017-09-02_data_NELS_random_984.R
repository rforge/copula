## By Marius Hofert

## Data preparation for the US National Education Longitudinal Study (NELS)
## of 1988 from Edward W. Frees


### 1 Cleaning #################################################################

## Read the (raw) data from Jed's website
raw <- read.csv("2017-09-02_data_NELS_random_984.csv")
str(raw)
stopifnot(length(unique(raw$schoolid)) == 489)

## Converting to the variables we want (are 'appended')
dat <- raw
dat <- within(dat, {
    ID       <- factor(schoolid)
    Math     <- math
    Science  <- sci
    Reading  <- read
    Minority <- factor(minority)
    SES      <- ses
    Female   <- factor(female)
    Public   <- factor(public)
    Size     <- schoolsize
    Urban    <- factor(urban)
    Rural    <- factor(rural)
    Filter   <- factor(filter_.) # => omitted later
})

## Grab out those we work with
dat <- dat[, c("ID", "Math", "Science", "Reading", "Minority", "SES", "Female",
               "Public", "Size", "Urban", "Rural")]
str(dat)
table(dat$Size) # => could be a factor, but we don't convert it

## Rural vs Urban
urban.rural <- cbind(dat$Urban, dat$Rural) - 1 # note: converted to 1:2 here!
table(rowSums(urban.rural))
## => At most one of 'urban' and 'rural' is 1, but they can both be 0


### 2 Save the cleaned data ####################################################

NELS88 <- dat
str(NELS88)
save(NELS88, file = "NELS88.rda", compress = "xz")

