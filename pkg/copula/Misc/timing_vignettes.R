## Idea, etc: By Marius Hofert

## Timing vignettes


setwd("../vignettes")

## Timing command
timer <- function(x) {
    tm <- system.time(res <- tryCatch(x, error = function(e) e))
    if(inherits(res, "error")) NA else round(tm[["elapsed"]], 2)
}
do.rm <- !(Sys.info()[["user"]] %in% c("maechler", "<you>"))#

## Determine all vignettes
all <- system("ls", intern = TRUE)
(all.Rnw <- all[grepl("\\.Rnw$", all)])
(all.Rmd <- all[grepl("\\.Rmd$", all)])
(n <- length(all.vign <- c(all.Rnw, all.Rmd)))
if(n < 1)
    stop("Found no vignettes in the current working directory")

## Time measurement
tm <- numeric(n)
for(j in seq_along(all.Rnw)) { # iterate over all .Rnw
    cat("==> File ", fil <- all.Rnw[j],": ")
    tm[j] <- timer(system(paste("R CMD Sweave --pdf", fil)))
    Rnwbasename <- tools::file_path_sans_ext(basename(fil))
    all <- system(paste0("ls ",Rnwbasename,"*"), intern = TRUE)
    if(do.rm) file.remove(all[!grepl(fil, all, fixed=TRUE)])# keep *.Rnw file (and backups of that)
    cat("\n \\--> finished ", fil, "\n-==-==-==-==-\n")
}
n1 <- length(all.Rnw)
## rmarkdown::render() :
## - Careful with environments: the code is executed in parent.frame(); see, ?render.
## - For example, if a global counter 'i' is used here, it is overwritten by wild_animals.Rmd.
## - If new.env() is used in render(), then dNAC.Rmd and AR_Clayton.Rmd fail with weird problems.
for(j in seq_along(all.Rmd)) { # iterate over all .Rmd
    cat("==> File ", fil <- all.Rmd[j],": ")
    tm[j+n1] <- timer(rmarkdown::render(fil))
    Rmdbasename <- tools::file_path_sans_ext(basename(fil))
    if(do.rm) file.remove(paste0(Rmdbasename,".html"))
    cat("\n-==-==-==-==-\n finished ", fil, "\n")
}
if(do.rm) file.remove("Rplots.pdf")

## Results
timings <- data.frame("Vignette" = all.vign,
                      "Elapsed.time.in.sec" = tm)
sink("../Misc/timings.out", split = TRUE)
timings[order(timings[,"Elapsed.time.in.sec"], decreasing = TRUE), ]
sink()
