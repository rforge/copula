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

## Check whether there are schools for which rural/urban differ
## and whether there are schools for which the size differs
num.uID <- table(dat$ID) # number of schools with a unique ID
uID <- as.numeric(names(num.uID)) # IDs of schools with unique ID
for(i in 1:length(uID)) {
    if(num.uID[i] >= 2) {
        dat. <- dat[dat$ID == uID[i],] # grab out all schools with that uID
        if((length(unique(dat.$Urban)) != 1) ||  # check whether Rural and Urban differ
           (length(unique(dat.$Rural)) != 1))
            cat("Urban or Rural of school with ID ",uID[i]," differ!\n")
        if(length(unique(dat.$Size)) != 1) # check whether the size of the school is the same for each unique ID
            cat("Provided sizes for school with ID ",uID[i]," differ!\n")
    }
}
## Result:
## - No such case => 'Urban' and 'Rural' could indicate the school or the student.
## - Also, 'Size' is indeed the same for each unique school (sanity check).


### 2 Save the cleaned data ####################################################

NELS88 <- dat
str(NELS88)
save(NELS88, file = "NELS88.rda", compress = "xz")

