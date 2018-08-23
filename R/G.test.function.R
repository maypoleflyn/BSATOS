

      G.test<-function (x, y = NULL, conservative = FALSE) 
{
    DNAME <- deparse(substitute(x))
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (is.matrix(x)) {
        if (min(dim(x)) == 1) 
            x <- as.vector(x)
    }
    if (!is.matrix(x) && !is.null(y)) {
        if (length(x) != length(y)) 
            stop("x and y must have the same length")
        DNAME <- paste(DNAME, "and", deparse(substitute(y)))
        OK <- complete.cases(x, y)
        x <- as.factor(x[OK])
        y <- as.factor(y[OK])
        if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
            stop("x and y must have at least 2 levels")
        x <- table(x, y)
    }
    if (any(x < 0) || any(is.na(x))) 
        stop("all entries of x must be nonnegative and finite")
    if ((n <- sum(x)) == 0) 
        stop("at least one entry of x must be positive")
    if (!is.matrix(x)) 
        stop("Could not make a 2-dimensional matrix")
    nrows <- nrow(x)
    ncols <- ncol(x)
    sr <- apply(x, 1, sum)
    sc <- apply(x, 2, sum)
    E <- outer(sr, sc, "*")/n
    g <- 0
    for (i in 1:nrows) {
        for (j in 1:ncols) {
            if (x[i, j] != 0) 
                g <- g + x[i, j] * log(x[i, j]/E[i, j])
        }
    }
    q <- 1
    if (conservative) {
        row.tot <- col.tot <- 0
        for (i in 1:nrows) {
            row.tot <- row.tot + 1/(sum(x[i, ]))
        }
        for (j in 1:ncols) {
            col.tot <- col.tot + 1/(sum(x[, j]))
        }
        q <- 1 + ((n * row.tot - 1) * (n * col.tot - 1))/(6 * 
            n * (ncols - 1) * (nrows - 1))
    }
    STATISTIC <- G <- 2 * g/q
    PARAMETER <- (nrow(x) - 1) * (ncol(x) - 1)
    PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
    if (!conservative) 
        METHOD <- "Log likelihood ratio (G-test) test of independence without correction"
    else METHOD <- "Log likelihood ratio (G-test) test of independence with Williams' correction"
    names(STATISTIC) <- "Log likelihood ratio statistic (G)"
    names(PARAMETER) <- "X-squared df"
    names(PVAL) <- "p.value"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD, data.name = DNAME, observed = x, 
        expected = E), class = "htest")
}

