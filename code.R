random.col <- function(r=runif(1), g=runif(1), b=runif(1), alpha=1)

## # Motivation section

## Suppose that our target distribution is a normal distribution with
## mean `m` and standard deviation `s`.  Obviously the mean of this
## distribution is `m`, but we can converge on that quickly.

## Parameters corresponding to a standard normal:
m <- 0
s <- 1

## Draw 10,000 samples from the distribution
set.seed(1)
samples <- rnorm(10000, m, s)

## The mean of the samples is very close to the true mean:
mean(samples)

## In fact, in this case, the expected error (the variance of the
## estimate) is 1/n, so we'd expect most values to lie within +/-
## 2/sqrt(n) = 0.02 of the true mean for 10,000 points.
summary(replicate(1000, mean(rnorm(10000, m, s))))

## This function computes the cumulative mean (that is, for element
## $k$, the sum of elements $1, 2, \ldots, k$ divided by $k$).
cummean <- function(x)
    cumsum(x) / seq_along(x)

## Here is the convergence towards the true mean (red line at 0).
plot(cummean(samples), type="l", xlab="Sample", ylab="Cumulative mean",
     panel.first=abline(h=0, col="red"), las=1)

## Transforming the x axis onto a log scale and showing another 30
## random approaches:
set.seed(1)
plot(cummean(samples), type="l", xlab="Sample", ylab="Cumulative mean",
     panel.first=abline(h=0, col="red"), las=1, log="x")
for (i in seq_len(30))
    lines(cummean(rnorm(10000, m, s)),
          col=rgb(runif(1), runif(1), runif(1), .5))

## True value:
a.true <- qnorm(p, m, s)

## Estimated by actual integration:
f <- function(x) dnorm(x, m, s)
g <- function(a)
    integrate(f, -Inf, a)$value
p <- 0.025
a.int <- uniroot(function(x) g(x) - p, c(-10, 0))$root

## Estimated by Monte Carlo integration -- has error associated with
## it!
a.mc <- unname(quantile(samples, p))
a.true - a.mc

a.mc <- replicate(100, unname(quantile(rnorm(10000, m, s), p)))
## Relative error on same order as the error around them mean.
summary(a.true - a.mc)

## # Markov chains

P <- rbind(c(.5,  .25, .25),
           c(.2,  .1,  .7),
           c(.25, .25, .5))

iterate.P <- function(x, P, n) {
    res <- matrix(NA, n+1, length(x))
    res[1,] <- x
    for (i in seq_len(n))
        res[i+1,] <- x <- x %*% P
    res
}

n <- 10
y1 <- iterate.P(c(1, 0, 0), P, n)
y2 <- iterate.P(c(0, 1, 0), P, n)
y3 <- iterate.P(c(0, 0, 1), P, n)

matplot(0:n, y1, type="l", lty=1, xlab="Step", ylab="y", las=1)
matlines(0:n, y2, lty=2)
matlines(0:n, y3, lty=3)

## which means that regardless of the starting distribution, there is
## a 32% chance of the chain being in state 1 after about 10 or more
## iterations *regardless of where it started*.

v <- eigen(t(P), FALSE)$vectors[,1]
v <- v/sum(v) # normalise eigenvector
## Points showing how close we are to convergence:
points(rep(10, 3), v, col=1:3)

## So, let's iterate the system, rather than the probability vector:
run <- function(i, P, n) {
    res <- integer(n)
    for (t in seq_len(n))
        res[[t]] <- i <- sample(nrow(P), 1, pr=P[i,])
    res
}

## Here's the chain running around...
samples <- run(1, P, 100)
plot(samples, type="s", xlab="Step", ylab="State", las=1)

plot(cummean(samples == 1), type="l", ylim=c(0, 1),
     xlab="Step", ylab="y", las=1)
lines(cummean(samples == 2), col=2)
lines(cummean(samples == 3), col=3)

## Run this out a little longer:
n <- 5000
set.seed(1)
samples <- run(1, P, n)
plot(cummean(samples == 1), type="l", ylim=c(0, 1), log="x",
     xlab="Step", ylab="y", las=1)
lines(cummean(samples == 2), col=2)
lines(cummean(samples == 3), col=3)
abline(h=v, lty=2, col=1:3)

## # MCMC sampling 1d

p <- 0.4
mu <- c(-1, 2)
sd <- c(.5, 2)
f <- function(x)
  p * dnorm(x, mu[1], sd[1]) + (1-p) * dnorm(x, mu[2], sd[2])
curve(f(x), col="red", -4, 8, n=301, las=1)

## Let's define a really simple minded proposal algorithm that samples
## from a normal distribution centred on the current point with a
## standard deviation of 4
q <- function(x) rnorm(1, x, 4)

## This implements the core algorithm:
step <- function(x, f, q) {
    ## Pick new point
    xp <- q(x)
    ## Plain Metropolis, so acceptance ratio is:
    alpha <- min(1, f(xp) / f(x))
    ## Accept new point with probability alpha:
    if (runif(1) < alpha)
        x <- xp
    ## Returning the point:
    x
}

## And this just takes care of running the MCMC for a number of
## steps:
run <- function(x, f, q, nsteps) {
    res <- matrix(NA, nsteps, length(x))
    for (i in seq_len(nsteps))
        res[i,] <- x <- step(x, f, q)
    drop(res)
}

## We pick a place to start (how about -10, just to pick a poor
## point).
res <- run(-10, f, q, 1000)
                
## Here are the first 1000 steps of the Markov chain, with the target
## density on the right:
layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
plot(res, type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1)
usr <- par("usr")
xx <- seq(usr[3], usr[4], length=301)
plot(f(xx), xx, type="l", yaxs="i", axes=FALSE, xlab="")

## Even with only a thousand (non-independent) samples, we're starting
## to resemble the target distribution fairly well.
hist(res, 50, freq=FALSE, main="", ylim=c(0, .4), las=1,
     xlab="x", ylab="Probability density")
z <- integrate(f, -Inf, Inf)$value
curve(f(x) / z, add=TRUE, col="red", n=200)

## Run for longer and things start looking a bunch better:
set.seed(1)
res.long <- run(-10, f, q, 50000)
hist(res.long, 100, freq=FALSE, main="", ylim=c(0, .4), las=1,
     xlab="x", ylab="Probability density", col="grey")
z <- integrate(f, -Inf, Inf)$value
curve(f(x) / z, add=TRUE, col="red", n=200)

## Now, run with different proposal mechanisms - one with a very wide
## standard deviation (33 units) and the other with a very small
## standard deviation (3 units).
res.fast <- run(-10, f, function(x) rnorm(1, x,  33), 1000)
res.slow <- run(-10, f, function(x) rnorm(1, x,  .3), 1000)

## Here is the same plot as above -- note the different ways that the
## three traces are moving around.
layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
plot(res, type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1,
     col="grey")
lines(res.fast, col="red")
lines(res.slow, col="blue")
plot(f(xx), xx, type="l", yaxs="i", axes=FALSE)

## The original (grey line) trace is bouncing around quite freely.

## In contrast, the red trace (large proposal moves) is suggesting
## terrible spaces in probability space and rejecting most of them.
## This means it tends to stay put for along time at once space.

## The blue trace proposes small moves that tend to be accepted, but
## it moves following a random walk for most of the trajectory.  It
## takes hundreds of iterations to even reach the bulk of the
## probability density.

## You can see the effect of different proposal steps in the
## autocorrelation among subsequent parameters -- these plots show the
## decay in autocorrelation coefficient between steps of different
## lags, with the blue lines indicating statistical independence.
par(mfrow=c(1, 3), mar=c(4, 2, 3.5, .5))
acf(res.slow, las=1, main="Small steps")
acf(res, las=1, main="Intermediate")
acf(res.fast, las=1, main="Large steps")

## From this, one can calculate the effective number of independent
## samples:
coda::effectiveSize(res)
coda::effectiveSize(res.fast)
coda::effectiveSize(res.slow)

## The chains both "mix" worse than that first one.

## This shows more clearly what happens as the chains are run for longer:
n <- 10^(2:5)
samples <- lapply(n, function(n) run(-10, f, q, n))
xlim <- range(sapply(samples, range))
br <- seq(xlim[1], xlim[2], length=100)

hh <- lapply(samples, function(x) hist(x, br, plot=FALSE))
ylim <- c(0, max(f(xx)))

## Showing 100, 1,000, 10,000 and 100,000 steps:
par(mfrow=c(2,2), mar=rep(.5, 4), oma=c(4, 4, 0, 0))
for (h in hh) {
    plot(h, main="", freq=FALSE, yaxt="n",
         ylim=range(h$density, ylim))
    curve(f(x), add=TRUE, col="red", n=300)
}

## # MCMC In two dimensions

## This is a function that makes a multivariate normal density given a
## vector of means (centre of the distribution) and
## variance-covariance matrix.
make.mvn <- function(mean, vcv) {
  logdet <- as.numeric(determinant(vcv, TRUE)$modulus)
  tmp <- length(mean) * log(2 * pi) + logdet
  vcv.i <- solve(vcv)

  function(x) {
    dx <- x - mean
    exp(-(tmp + rowSums((dx %*% vcv.i) * dx))/2)
  }
}

## To make it interesting, let's define a function that is the sum of
## two mvns:

mu1 <- c(-1, 1)
mu2 <- c(2, -2)
vcv1 <- matrix(c(1, .25, .25, 1.5), 2, 2)
vcv2 <- matrix(c(2, -.5, -.5, 2), 2, 2)
f1 <- make.mvn(mu1, vcv1)
f2 <- make.mvn(mu2, vcv2)
f <- function(x)
    f1(x) + f2(x)

x <- seq(-5, 6, length=71)
y <- seq(-7, 6, length=61)
xy <- expand.grid(x=x, y=y)
z <- matrix(apply(as.matrix(xy), 1, f), length(x), length(y))

image(x, y, z, las=1)
contour(x, y, z, add=TRUE)

## Sampling from multivariate normals is also fairly straightforward,
## but we'll draw samples from this using mcmc.

## There are a bunch of different strategies here -- we could propose
## moves in both dimensions simultaneously, or we could sample along
## each axis independently.  Both strategies will work, though they
## will again differ in how rapidly they mix.

## Assume that we don't actually know how to sample from a mvn (it's
## not actually hard, but this is simpler), let's make a proposal
## distribution that is uniform in two dimensions, sampling from the
## square with width 'd' on each side.
q <- function(x, d=8)
    x + runif(length(x), -d/2, d/2)

x0 <- c(-4, -4)
set.seed(1)
samples <- run(x0, f, q, 1000)

image(x, y, z, xlim=range(x, samples[,1]), ylim=range(x, samples[,2]))
contour(x, y, z, add=TRUE)
lines(samples[,1], samples[,2], col="#00000088")

## Drawing a ton of samples"
set.seed(1)
samples <- run(x0, f, q, 100000)

## Compare the sampled distribution against the known distribution:
smoothScatter(samples)
contour(x, y, z, add=TRUE)

### library(hexbin)
### hbin <- hexbin(samples[,1], samples[,2], xbins = 40)
### plot(hbin, legend=FALSE)
### contour(x, y, z, add=TRUE)

## Then we can do things that are really nice to do here, but hard to
## to in general.  For example, what is the *marginal distribution* of
## parameter 1:
hist(samples[,1], freq=FALSE, main="", xlab="x",
     ylab="Probability density")

## Computing this properly is tricky - we need to integrate over all
## possible values of the second parameter for each value of the
## first.  Then, because the target function is not itself normalised,
## we have to divide that through by the value of integrating over the
## first dimension (this is the total area under the distribution).
m <- function(x1) {
    g <- Vectorize(function(x2) f(c(x1, x2)))
    integrate(g, -Inf, Inf)$value
}

xx <- seq(min(samples[,1]), max(samples[,1]), length=201)
yy <- sapply(xx, m)
z <- integrate(splinefun(xx, yy), min(xx), max(xx))$value

hist(samples[,1], freq=FALSE, main="", las=1, xlab="x", 
     ylab="Probability density", ylim=c(0, 0.25))
lines(xx, yy/z, col="red")
