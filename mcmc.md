\newcommand{\ud}{\mathrm{d}}

# What is MCMC and when would you use it?

## MCMC is a way of sampling from a distribution

It's only one of many algorithms for doing so.  The term stands for
"Markov Chain Monte Carlo", because it is a type of "Monte Carlo"
(i.e., a random) method that uses "Markov chains" (we'll discuss these
later).  MCMC is just one type of Monte Carlo method, and we'll
describe these first.

## Why would I want to sample from a distribution?

You may not realise you want to (and really, you may not actually want
to).  However, sampling from a distribution turns out to be the
easiest way of solving some problems.  

Probably the most common way that MCMC is used is to draw samples from
the **posterior distribution** of some model in Bayesian inference.
With these samples, you can then ask things like "what is the mean and
credibility interval for a parameter".

For example, suppose that you have fit a model where the posterior
probability density is some function $f$ of parameters $x, y$.  Then,
to compute the mean value of parameter $x$, you would compute

$$
\bar x = \int\int x f(x, y)~\ud x~\ud y
$$

which you can read simply as "the value of $x$ multiplied by the
probability of parmeters $(x, y)$, integrated over all possible values
that $x$ and $y$ could take.

An alternative way to compute this value is simulate $k$ observations:
$\{(x,y)^{(1)}, \ldots, (x,y)^{(k)}\}$ from $f(x,y)$ and compute the sample
mean:

$$\bar x \approx \frac{1}{k} \sum_j x^{(j)}$$

where $x^{(j)}$ is the the $x$ from the $j$th sample from the
distribution.

If these samples are independent samples from the distribution, then
as $k \to \infty$ the estimated mean of $x$ will converge on the true
mean.

*See `code.R`, Motivation section for code to use here*

How does this work?  Consider the integral

$$
\int_a^b h(x) \ud x
$$

If this can be decomposed into the product of a function $f(x)$ and a
probability density function $p(x)$, then

$$
\int_a^b h(x) \ud x = \int_a^b f(x)p(x) \ud x
$$

Note that the right hand side is simply the expectation $E[f(x)]$.  By
the Law of Large Numbers", the expected value is the limit of the sample
mean as the sample size grows to infinity".  So we can approximate
$E[f(x)]$ as

$$
\frac{1}{n}\sum_{i=1}^n f(x_i).
$$

You can do lots of similar things with these.  For example, if you
want to draw a 95\% credibility interval around the estimate $\bar x$,
you could estimate the bottom component of that by solving

$$
0.025 = \int_{-\infty}^a\int x f(x, y)~\ud y~\ud x
$$

for $a$.  Or, you can just take the sample quantile from your series
of sampled points.

*See second part of `code.R` here*

## Why doesn't "normal statistics" use Monte Carlo methods?

For many problems in traditionally taught statistics, rather than
sampling from a distribution you **maximise or maximise a function**.
So we'd take some function that describes the likelihood and maximise
it (maximum likelihood inference), or some function that computes the
sum of squares and minimise it.


The reasons for this difference are a little subtle, but boil down to
whether or not you feel that you could possibly put a probability
distribution over a parameter -- is it something you could sample?
Fisher in particular had strong thoughts on this, thoughts which are
argued more recently by AWF Edwards in the book likelihood.  To avoid
having to sample from a distribution (or really, to avoid the idea
that one could draw samples from a probability distribution of
parameters), error estimates in frequentist statistics tend to either
be asymptotic large-data estimates or perhaps bootstrap based.

However, every time you are doing something like a "least squares"
analysis you are minimising a function -- in this case the sum of
squared differences between your predicted model and your data.  The
"black box" that forms the optimisation proceedure is the analagous
step do doing an MCMC.

# Markov Chain Monte Carlo

At this point, suppose that there is some target distribution that
we'd like to sample from, but that we cannot just draw independent
samples from like we did before.  There is a solution for doing this
using the Markov Chain Monte Carlo (MCMC).  First, we have to define
some things so that the next sentence makes sense: *What we're going
to do is try to construct a Markov chain that has our
hard-to-sample-from target distribution as its stationary
distribution*.  This section is going to skip all the gory details,
and just show the mechanics.  We'll cover how it works later.

## Definitions

Let $X_t$ denote the value of some random variable at time $t$.  A
Markov chain generates a series of samples $\{X_0, X_1, ..., X_t\}$.

Markov chains satisfy the *Markov property*.  The Markov property is
the probabilistic version of "what happens in Vegas stays in Vegas";
basically it doesn't matter how you got to some state $x$, the
probability of transition out of $x$ is unchanged, or:

$$\Pr(X_{t+1} = x | X_t = x_t, X_{t-1} = x_{t-1}, \ldots, X_0 = x_0) =
\Pr(X_{t+1} = x | X_t = x_t).$$

The transition from one step to the next is described by the
*transition kernel*, which can be described by the probability (or for
continuous variables the probability *density*) of a transition from
state $i$ to state $j$ as 

$$P(i \to j) = \Pr(X_{t+1} = x_j | X_t = x_i).$$

If the process is *irreducible* (every state is visitable from every
other state) and *aperiodic* (the number of steps between two visits
of a state is not a fixed integer multiple number of steps) then it
has a *stationary distribution*.

Let $pi_j(t) = \Pr(X_t = s_j)$ be the probability that the chain is in
state $j$ at time (step) $t$, and define $\vec\pi(t)$ be the vector of
probabilites over possible states.  Then, given $\vec\pi(t)$, we can
compute $\vec\pi(t+1)$ using the *Chapman-Kolmogorov* equation.

$$\pi_i(t+1) = \sum_k \pi_k(t) \Pr(k \to i),$$

that is; the probability that we were in state $k$ multiplied by the
probability of making the transition from $k$ to $i$, summed over all
possible source states $k$.  Using the book-keeping of linear algebra,
let $\mathbf{P}$ be the *probability transition matrix* -- the matrix
whose $i,j$th element is $P(i \to j)$, and rewrite the above equation
as

$$\vec\pi(t + 1) = \vec\pi(t)\mathbf{P}$$

Note that we can iterate this equation easily:

$$\vec\pi(t+2) = \vec\pi(t+1)\mathbf{P}$$
$$\vec\pi(t+2) = \vec\pi(t)\mathbf{P}\mathbf{P}$$
$$\vec\pi(t+2) = \vec\pi(t)\mathbf{P}^2$$

If there is some vector $\vec\pi^*$ that satisfies

$$\vec\pi^* = \vec\pi^*\mathbf{P}$$

then $\vec\pi^*$ is the *stationary distribution* of this Markov
chain.  Mathematically, $\vec\pi^*$ is the left eigenvector assicated
with the eigenvalue = 1.

*See `code.R` here, section Markov chains*

A sufficient (but not necessary) condition for the existance of a
stationary distribution is Detailed Balance, which says:

$$P(j \to k) \pi_j^* = P(k \to j) \pi_k^*$$

This imples that the chain is *reversible*.  The reason why this
condition implies that a stationary distribution exists is that it
implies

$$\vec\pi^* = \vec\pi^*\mathbf{P}$$

Summing both sides of the detailed balance equation over states $j$

$$\sum_j P(j \to k) \pi_j^* = \sum_j P(k \to j) \pi_k^*$$

The term on the left is equal to the $k$th element of
$\vec\pi^*\mathbf{P}$ and the term on the right can be factored:

$$\sum_j P(k \to j) \pi_k^* = \pi_k^* \sum_j P(k \to j)$$

Then, because $\sum_j P(k \to j) = 1$ (because $P$ is a transition
probability function, by the law of total probability things go
*somewhere* with probability 1), so the right hand side is $\pi_k^*$,
so we have 

$$(\vec\pi^*\mathbf{P})_k = \pi_k^*$$

which holds for all $k$ so

$$\vec\pi^*\mathbf{P} = \vec\pi^*$$

## The Metropolis algorithm

This is the simplest MCMC algorithm.  This section is not intended to
show how to design efficient MCMC samplers, but just to see that they
do in fact work.  What we're going to do is have some distribution
that we want to sample from, and we're going to be able to evaluate
some function $f(x)$ that is *proportional* to the probability density
of the target distribution (that is, if $p(x)$ is the probability
density function itself, $f(x) \propto p(x)$, i.e., $f(x) = p(x) / Z$,
where $Z = \int f(x) \ud x$).  Note that $x$ might be a vector or a
scalar.

We also need a probability density function $P$ that we *can* draw
samples from.  For the simplest algorithm, this **proposal**
distribution is symmetric, that is $P(x\to x^\prime) = P(x^\prime \to
x)$.

The algorithm proceeds as follows.

1. Start in some state $x_t$.
2. Propose a new state $x^\prime$
3. Compute an acceptance probability 
$$\alpha = \min\left[1, \frac{f(x^\prime)}{f(x)}\right]$$
4. Draw some uniformly distributed random number $u$ from $[0,1]$; if
$u < \alpha$ accept the point, setting $x_{t+1} = x^\prime$.
Otherwise reject it and set $x_{t+1} = x_t$.

Note that in step 3. above, the unknown normalising constant drops out
because

$$\frac{p(x^\prime)}{p(x)} = \frac{f(x^\prime)}{Z} \frac{Z}{f(x)} =
\frac{f(x^\prime)}{f(x)}$$

This will generate a series of samples $\{x_0, x_1, \ldots\}$.  Note
that where the proposed sample is rejected, the same *value* will be
present in consecutive samples.

*See `code.R`, MCMC sampling 1d here*

## MCMC sampling in 1d (single parameter) problems

Here is a target distribution to sample from.  It's the weighted sum
of two normal distributions.  This sort of distribution is fairly
straightforward to sample from, but let's draw samples with MCMC.  The
probability density function is

$$
f(x) = \frac{1}{z}
(w_1 \exp(x - \mu_1)^2 / (2 \sigma_1^2) +
w_2 \exp(x - \mu_2)^2 / (2 \sigma_2^2))
$$
