

library(stats4) ## for mle
library(chipseq)
library(lattice)


rawData <- function(data, chr, start, end)
{
    data <- data[[chr]]
    data <- lapply(data, function(x) x[x >= start & x <= end])
    data
}

## sigma = 0 version 

loglik <- function(theta, p, mu, sigma = 0, raw.data, width)
{
    if (p < 0 || p > 1 || mu < 0) return(-Inf)
    f0 <- function(x) 1/width
    f1.minus <- function(x) ifelse(x >= theta & x <= theta + mu, 1/mu, 0)
    f1.plus <- function(x) ifelse(x <= theta & x >= theta - mu, 1/mu, 0)
    sum(log(p * (f0(raw.data[["+"]])) + (1-p) * f1.plus(raw.data[["+"]]))) + 
        sum(log(p * (f0(raw.data[["-"]])) + (1-p) * f1.minus(raw.data[["-"]])))
}

negllik.p.mu <- function(theta, sigma = 0, raw.data, width)
{
    function(p, mu)
        -loglik(theta = theta, p = p, mu = mu, sigma = sigma, raw.data = raw.data, width = width)
}

doTheta <- function(theta, raw.data, width)
{
    lfun <- negllik.p.mu(theta = theta, raw.data = raw.data, width = width)
    m <- mle(lfun, start = list(p = 0.5, mu = 200))
    c(M2LL = 2 * m@min, m@coef, sd = sqrt(diag(m@vcov)))
}

doRegion <- function(data, chr, start, end, plot = interactive(), by = 5)
{
    rd <- rawData(data = data, chr = chr, start = start, end = end)
    w <- end - start
    doTheta <- function(theta)
    {
        lfun <- negllik.p.mu(theta = theta, raw.data = rd, width = w)
        m <- try(mle(lfun, start = list(p = 0.5, mu = 200)), silent = TRUE)
        if (inherits(m, "try-error")) rep(NA_real_, 5) 
        else c(M2LL = 2 * m@min, m@coef, sd = sqrt(diag(m@vcov)))
    }
    i <- seq(start, end, by = by)
    ans1 <- 
        sapply(i, 
               function(theta) {
                   if (interactive())
                       cat("\r", theta, "/",
                           end-theta,
                           "         ")
                   doTheta(theta)
               })
    if (interactive()) cat("\n")
    ans1 <- as.data.frame(t(ans1))
    names(ans1) <- c("M2LL", "p", "mu", "sd.p", "sd.mu")
    ans1$theta <- i
    ## do region close to min with finer resolution
    argmin.theta <- ans1$theta[which.min(ans1$M2LL)]
    extend <- 20
    i <- seq(argmin.theta - extend, argmin.theta + extend, by = 1)
    ans2 <- 
        sapply(i, 
               function(theta) {
                   if (interactive())
                       cat("\r", theta, "/",
                           argmin.theta + extend -theta,
                           "         ")
                   doTheta(theta)
               })
    if (interactive()) cat("\n")
    ans2 <- as.data.frame(t(ans2))
    names(ans2) <- c("M2LL", "p", "mu", "sd.p", "sd.mu")
    ans2$theta <- i
    ans <- rbind(subset(ans1, theta < min(i)),
                 ans2,
                 subset(ans1, theta > max(i)))
    argmin.theta <- ans$theta[which.min(ans$M2LL)]
    ## ans$plus <- density(rd[["+"]], from = start, to = end, n = length(i))$y
    ## ans$minus <- density(rd[["-"]], from = start, to = end, n = length(i))$y
    if (plot)
    {
        plot(stripplot(which ~ data, do.call(make.groups, rd),
                       panel = function(...) {
                           panel.abline(v = argmin.theta)
                           panel.stripplot(...)
                       },
                       jitter = TRUE, xlab = "",
                       main = sprintf("%s [ %g - %g ]", chr, start, end),
                       scales = list(x = list(draw = FALSE)),
                       xlim = range(ans$theta)),
             position = c(0, 0.7, 1, 1))
        plot(xyplot(M2LL + p + mu ~ theta, data = ans, type = "l", outer = TRUE,
                    panel = function(...) {
                        panel.abline(v = argmin.theta)
                        panel.xyplot(...)
                    },
                    lwd = 3,
                    scales = list(y = list(relation = "free", rot = 0)), strip = FALSE, strip.left = TRUE, 
                    ylab = "", layout = c(1, 3), xlim = range(ans$theta)),
             position = c(0, 0, 1, 0.7), newpage = FALSE)
    }
    ## invisible(ans)
    ans[which.min(ans$M2LL), ]
}


load("myodMyo.rda")

ctubes <- combineLaneReads(myodMyo[c("2","4","7")])

peakIslands <- function(chr, data, g = extendReads(data[[chr]]), lower = 10)
{
    if (interactive()) message(chr)
    s <- slice(coverage(g, 1, max(end(g))), lower = 1)
    s <- s[viewMaxs(s) >= lower]
    data.frame(chrom = chr, start = start(s), end = end(s),
               maxs = viewMaxs(s), sums = viewSums(s),
               stringsAsFactors = FALSE)
}

peaks <-
    sapply(names(ctubes),
           peakIslands, data = ctubes, lower = 12,
           simplify = FALSE)


doRow <- function(peakSet, i = 1, data, ...)
{
    x <- peakSet[i,]
    doRegion(data = data, chr = x$chrom, 
             start = x$start, end = x$end,
             ...)
}

pdf("chr14_peak_mle.pdf", width = 7, height = 9)

chr14.summary <- 
    lapply(seq_len(nrow(peaks$chr14)),
           function(i) {
               message(i, "/", nrow(peaks$chr14))
               doRow(peaks$chr14, i, data = ctubes, plot = TRUE)
           })

dev.off()


chr14.summary <- do.call(rbind, chr14.summary)
save(chr14.summary, file = "chr14.summary.rda")


xyplot(p ~ mu, chr14.summary)


