## By Marius Hofert

## Timing vignettes


setwd("../vignettes")

## Timing command
timer <- function(x) {
    tm <- system.time(res <- tryCatch(x, error = function(e) e))
    if(is(res, "simpleError")) NA else round(tm[["elapsed"]], 2)
}

## Determine all vignettes
all <- system("ls", intern = TRUE)
all.Rnw <- all[grepl(".Rnw", all)]
all.Rmd <- all[grepl(".Rmd", all)]
all.vign <- c(all.Rnw, all.Rmd)
n <- length(all.vign)
if(n < 1)
    stop("Found no vignettes in the current working directory")

## Time measurement
## Careful with environments: the code is executed in parent.frame(); see, e.g., ?render.
## For example, if a global counter 'i' is used here, it is overwritten by wild_animals.Rmd.
## If new.env() is used in render(), then dNAC.Rmd and AR_Clayton.Rmd fail with weird problems.
tm <- c()
for(j in seq_along(all.Rnw)) { # iterate over all .Rnw
    tm <- c(tm, timer(system(paste("R CMD Sweave --pdf", all.Rnw[j]))))
    Rnwbasename <- tools::file_path_sans_ext(basename(all.Rnw[j]))
    all <- system(paste0("ls ",Rnwbasename,"*"), intern = TRUE)
    all <- all[all != all.Rnw[j]] # keep original file
    file.remove(all)
}
for(j in seq_along(all.Rmd)) { # iterate over all .Rmd
    tm <- c(tm, timer(rmarkdown::render(all.Rmd[j])))
    Rmdbasename <- tools::file_path_sans_ext(basename(all.Rmd[j]))
    file.remove(paste0(Rmdbasename,".html"))
}
file.remove("Rplots.pdf")

## Results
timings <- data.frame("Vignette" = all.vign, "Elapsed.time.in.sec" = tm)
timings[order(timings[,"Elapsed.time.in.sec"], decreasing = TRUE), ]
