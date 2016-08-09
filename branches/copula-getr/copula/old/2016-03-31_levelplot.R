## This functionality is now provided by contourplot2()

##' @title A Level Plot with Nice Defaults
##' @param x A numeric matrix or as.matrix(.)able
##' @param aspect The aspect ratio
##' @param xlim The x-axis limits
##' @param ylim The y-axis limits
##' @param xlab The x-axis label
##' @param ylab The y-axis label
##' @param cuts The number of levels
##' @param scales See ?levelplot
##' @param par.settings Additional arguments passed to 'par.settings' (some are set)
##' @param contour A logical indicating whether contour lines should be plotted
##' @param ... Further arguments passed to levelplot()
##' @return A levelplot() object
##' @author Marius Hofert
levelplot2 <- function(x, aspect = 1,
                       xlim = extendrange(x[,1], f = 0.04),
                       ylim = extendrange(x[,2], f = 0.04),
                       xlab = NULL, ylab = NULL,
                       cuts = 20, scales = list(alternating = c(1,1), tck = c(1,0)),
                       par.settings = list(regions = list(col = gray(seq(0.4, 1, length.out=cuts+1)))),
                       contour = TRUE, ...)
{
    ## Checking
    if(!is.matrix(x)) x <- as.matrix(x)
    if(ncol(x) != 3) stop("'x' should be trivariate")

    ## Labels
    if(is.null(xlab) || is.null(ylab)) {
        colnms <- colnames(x)
        if(sum(nzchar(colnms)) != 2) {
            xlab <- expression(u[1])
            ylab <- expression(u[2])
        } else { # 'x' has column names => parse them
            xlab <- parse(text = colnms[1])
            ylab <- parse(text = colnms[2])
        }
    }

    ## Level plot
    levelplot(x[,3] ~ x[,1] * x[,2], aspect = aspect, xlim = xlim, ylim = ylim,
              xlab = xlab, ylab = ylab, cuts = cuts,
              scales = scales, par.settings = par.settings,
              contour = contour, ...)
}

