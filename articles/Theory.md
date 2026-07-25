# Theory

## Theory the `gips` is based on

The package is based on the article
[\[1\]](https://arxiv.org/abs/2004.03503). There the math behind the
package is precisely demonstrated, and all the theorems are proven.

In this vignette, we would like to give a gentle introduction. We want
to point out all the most important results from this work from the
user’s point of view. We will also show examples of those results in the
`gips` package.

``` r

library(gips)
```

As mentioned in the abstract, the outline of the paper is to “derive the
distribution of the maximum likelihood estimate of the covariance
parameter $`\Sigma`$ (…)” and then to “perform Bayesian model selection
in the class of complete Gaussian models invariant by the action of a
subgroup of the symmetric group (…)”. Those ideas are implemented in the
`gips` package.

## Basic definitions

Let $`V=\{1,\ldots,p\}`$ be a finite index set, and for every
$`i\in V`$, $`Z^{(i)}=(Z_1^{(i)},\ldots, Z_p^{(i)})^\top`$ be a
multivariate random variable following a centered Gaussian model
$`\mathrm{N}_p(0,\Sigma)`$, and let $`Z^{(1)},\ldots, Z^{(n)}`$ be an
i.i.d. (independent and identically distributed) sample from this
distribution. Name the whole sample $`Z = (Z^{(1)},\ldots, Z^{(n)})`$.

Let $`\mathfrak{S}_p`$ denote the symmetric group on $`V`$, that is, the
set of all permutations on $`\{1,\ldots,p\}`$ with function composition
as the group operation. Let $`\Gamma`$ be an arbitrary subgroup of
$`\mathfrak{S}_p`$. The model $`\mathrm{N}_p(0,\Sigma)`$ is said to be
invariant under the action of $`\Gamma`$ if for all $`g\in \Gamma`$,
$`g\cdot\Sigma\cdot g^\top=\Sigma`$ (here, we identify a permutation
$`g`$ with its permutation matrix).

For a subgroup $`\Gamma \subset  \mathfrak{S}_p`$, we define the colored
space, i.e. the space of symmetric matrices invariant under $`\Gamma`$,
``` math
\mathcal{Z}_{\Gamma} := \{S \in \mathrm{Sym}(p;\mathbb{R})\colon S_{i,j} = S_{\sigma(i),\sigma(j)} \text{ for all }\sigma \in \Gamma\mbox{ for all }i,j\in V\},
```
and the colored cone of positive definite matrices valued in
$`\mathcal{Z}_{\Gamma}`$,
``` math
\mathcal{P}_{\Gamma} := \mathcal{Z}_{\Gamma} \cap \mathrm{Sym}^+(p;\mathbb{R}).
```

## Block Decomposition - \[1\], Theorem 1

The main theoretical result in this theory (Theorem 1 in
[\[1\]](https://arxiv.org/abs/2004.03503)) states that given a
permutation subgroup $`\Gamma`$ there exists an orthogonal matrix
$`U_\Gamma`$ such that all the symmetric matrices
$`S\in\mathcal{Z}_\Gamma`$ can be transformed into block-diagonal form.

The exact form of blocks depends on so-called *structure constants*
$`(k_i,d_i,r_i)_{i=1}^L`$. It is worth pointing out that constants
$`k = d`$ for cyclic group $`\Gamma = \left<\sigma\right>`$ and that
`gips` searches within cyclic subgroups only.

#### Examples

``` r

p <- 6
S <- matrix(c(1.1,0.9,0.8,0.7,0.8,0.9,
              0.9,1.1,0.9,0.8,0.7,0.8,
              0.8,0.9,1.1,0.9,0.8,0.7,
              0.7,0.8,0.9,1.1,0.9,0.8,
              0.8,0.7,0.8,0.9,1.1,0.9,
              0.9,0.8,0.7,0.8,0.9,1.1), nrow = p)
```

![](Theory_files/figure-html/th1_2-1.png)

`S` is a symmetric matrix invariant under the group
$`\Gamma = \left<(1,2,3,4,5,6)\right>`$.

``` r

g_perm <- gips_perm("(1,2,3,4,5,6)", p)
U_Gamma <- prepare_orthogonal_matrix(g_perm)

block_decomposition <- t(U_Gamma) %*% S %*% U_Gamma
round(block_decomposition, 5)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  5.2  0.0  0.0  0.0  0.0  0.0
#> [2,]  0.0  0.5  0.0  0.0  0.0  0.0
#> [3,]  0.0  0.0  0.5  0.0  0.0  0.0
#> [4,]  0.0  0.0  0.0  0.1  0.0  0.0
#> [5,]  0.0  0.0  0.0  0.0  0.1  0.0
#> [6,]  0.0  0.0  0.0  0.0  0.0  0.2
```

![](Theory_files/figure-html/th1_4-1.png)

The transformed matrix is in the block-diagonal form of
[\[1\]](https://arxiv.org/abs/2004.03503), Theorem 1. Blank entries are
off-block entries and equal to 0. The result was rounded to the 5th
place after the decimal to hide the inaccuracies of floating point
arithmetic.

Let’s see the other example:

``` r

p <- 6
S <- matrix(c(1.2,0.9,0.9,0.4,0.2,0.1,
              0.9,1.2,0.9,0.1,0.4,0.2,
              0.9,0.9,1.2,0.2,0.1,0.4,
              0.4,0.1,0.2,1.2,0.9,0.9,
              0.2,0.4,0.1,0.9,1.2,0.9,
              0.1,0.2,0.4,0.9,0.9,1.2), nrow = p)
```

![](Theory_files/figure-html/th1_6-1.png)

Now, `S` is a symmetric matrix invariant under the group
$`\Gamma = \left<(1,2,3)(4,5,6)\right>`$.

``` r

g_perm <- gips_perm("(1,2,3)(4,5,6)", p)
U_Gamma <- prepare_orthogonal_matrix(g_perm)

block_decomposition <- t(U_Gamma) %*% S %*% U_Gamma
round(block_decomposition, 5)
#>      [,1] [,2]   [,3]    [,4]    [,5]   [,6]
#> [1,]  3.0  0.7 0.0000  0.0000  0.0000 0.0000
#> [2,]  0.7  3.0 0.0000  0.0000  0.0000 0.0000
#> [3,]  0.0  0.0 0.3000  0.0000  0.2500 0.0866
#> [4,]  0.0  0.0 0.0000  0.3000 -0.0866 0.2500
#> [5,]  0.0  0.0 0.2500 -0.0866  0.3000 0.0000
#> [6,]  0.0  0.0 0.0866  0.2500  0.0000 0.3000
```

![](Theory_files/figure-html/th1_8-1.png)

Again, this result is in accordance with
[\[1\]](https://arxiv.org/abs/2004.03503), Theorem 1. Notice the zeros
in `block_decomposition`:
``` math
\forall_{i\in\{1,2\},j\in\{3,4,5,6\}}\text{block_decomposition}[i,j] = 0
```

## Project Matrix - \[1, Eq. (6)\]

One can also take any symmetric square matrix `S` and find the
orthogonal projection on $`\mathcal{Z}_{\Gamma}`$, the space of matrices
invariant under the given permutation:

``` math
\pi_\Gamma(S) := \frac{1}{|\Gamma|}\sum_{\sigma\in\Gamma}\sigma\cdot S\cdot\sigma^\top
```

The projected matrix is the element of the cone
$`\pi_\Gamma(S)\in\mathcal{Z}_{\Gamma}`$, which means:
``` math
\forall_{i,j\in \{1,\ \dots,\ p\}} \pi_\Gamma(S)[i,j] = \pi_\Gamma(S)[\sigma(i),\sigma(j)] \text{ for all }\sigma\in\Gamma
```

So it has some identical elements.

#### Trivial case

Note that for $`\Gamma = \{\text{id}\} = \{(1)(2)\dots(p)\}`$ we have
$`\pi_{\{\text{id}\}}(S) = S`$.

So no additional assumptions are made; thus the standard covariance
estimator is the best we can do.

#### Example

Let `S` be any symmetric square matrix:

``` r

round(S, 2)
#>        [,1]  [,2]   [,3]  [,4]   [,5]  [,6]
#> [1,] 108.99  3.44 -14.41 14.90 -15.28 -5.78
#> [2,]   3.44 27.71  -3.45  2.56   4.81 -4.46
#> [3,] -14.41 -3.45  42.59  8.15  -7.57 -0.95
#> [4,]  14.90  2.56   8.15 29.94  -1.52  0.83
#> [5,] -15.28  4.81  -7.57 -1.52  10.85  5.56
#> [6,]  -5.78 -4.46  -0.95  0.83   5.56 27.04
```

![](Theory_files/figure-html/def3_2-1.png)

One can project this matrix, for example, on
$`\Gamma = \left<(1,2)(3,4,5,6)\right>`$:

``` r

S_projected <- project_matrix(S, perm = "(1,2)(3,4,5,6)")
round(S_projected, 2)
#>       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
#> [1,] 68.35  3.44 -7.90  2.62 -7.90  2.62
#> [2,]  3.44 68.35  2.62 -7.90  2.62 -7.90
#> [3,] -7.90  2.62 27.60  2.81 -3.37  2.81
#> [4,]  2.62 -7.90  2.81 27.60  2.81 -3.37
#> [5,] -7.90  2.62 -3.37  2.81 27.60  2.81
#> [6,]  2.62 -7.90  2.81 -3.37  2.81 27.60
```

![](Theory_files/figure-html/def3_4-1.png)

Notice in the `S_projected` matrix, there are identical elements
according to the equation from the beginning of this section. For
example, `S_projected[1,1] = S_projected[2,2]`.

### $`C_\sigma`$ and `n0`

It is a well-known fact that without additional assumptions, the Maximum
Likelihood Estimator (MLE) of the covariance matrix in the Gaussian
model exists if and only if $`n \ge p`$. However, if the additional
assumption is added as the covariance matrix is invariant under
permutation $`\sigma`$, then the sample size $`n`$ required for the MLE
to exist is lower than $`p`$. It is equal to the number of cycles,
denoted hereafter by $`C_\sigma`$.

For example, if the permutation $`\sigma = (1,2,3,4,5,6)`$ is discovered
by the
[`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md)
function, then there is a single cycle in it $`C_\sigma = 1`$. Therefore
a single observation would be enough to estimate a covariance matrix
with
[`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md).
If the permutation $`\sigma = (1,2)(3,4,5,6)`$ is discovered, then
$`C_\sigma = 2`$ and so 2 observations would be enough.

To get this $`C_\sigma`$ number in `gips`, one can call
[`summary()`](https://rdrr.io/r/base/summary.html) on the appropriate
`gips` object:

``` r

g1 <- gips(S, n, perm = "(1,2,3,4,5,6)", was_mean_estimated = FALSE)
summary(g1)$n0
#> [1] 1
g2 <- gips(S, n, perm = "(1,2)(3,4,5,6)", was_mean_estimated = FALSE)
summary(g2)$n0
#> [1] 2
```

This is called `n0` and not $`C_\sigma`$, because it is increased by 1
when the mean was estimated:

``` r

S <- cov(Z)
g1 <- gips(S, n, perm = "(1,2,3,4,5,6)", was_mean_estimated = TRUE)
summary(g1)$n0
#> [1] 2
g2 <- gips(S, n, perm = "(1,2)(3,4,5,6)", was_mean_estimated = TRUE)
summary(g2)$n0
#> [1] 3
```

## Bayesian model selection

When one has the data matrix `Z`, one would like to know if it has a
hidden structure of dependencies between features. Luckily, the paper
demonstrates a way how to find it.

#### General workflow

1.  Choose the prior distribution on $`\Gamma`$ and $`\Sigma`$.
2.  Calculate the posteriori distribution (up to a normalizing constant)
    by formula [\[1\]](https://arxiv.org/abs/2004.03503), (30).
3.  Use the Metropolis-Hastings algorithm to find the permutation with
    the biggest value of the posterior probability
    $`\mathbb{P}(\Gamma|Z)`$.

#### Details on the prior distribution

The considered prior distribution of $`\Gamma`$ and $`K=\Sigma^{-1}`$:

1.  $`\Gamma`$ is uniformly distributed on the set of all cyclic
    subgroups of $`\mathfrak{S}_p`$.
2.  $`K`$ given $`\Gamma`$ follows the Diaconis-Ylvisaker conjugate
    prior distribution with parameters $`\delta`$ (real number,
    $`\delta > 2`$) and $`D`$ (symmetric, positive definite square
    matrix of the same size as `S`), see
    [\[1\]](https://arxiv.org/abs/2004.03503), Sec. 3.4.

Footnote: Actually in some cases smaller (but still positive) values for
$`\delta`$ parameter are theoretically correct. Refer to the
[\[1\]](https://arxiv.org/abs/2004.03503).

#### `gips` technical details

In `gips`, $`\delta`$ is named `delta`, and $`D`$ is named `D_matrix`.
By default, they are set to $`3`$ and `diag(p)`, respectively. However,
it is worth running the procedure for several parameters `D_matrix` of
form $`C\cdot diag(p)`$ for positive constant $`C`$. Large $`C`$
(compared to the data) favors small groups.

One can calculate the logarithm of formula (30) with function
[`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).

#### Interpretation

When all assumptions are met, the formula (30) puts a number on each
permutation’s cyclic group. The bigger its value, the more likely the
data was drawn from that model.

When one finds the permutations group $`c_{\text{max}}`$ that maximizes
(30),
``` math
c_{\text{map}} = \operatorname{arg\,max}_{c\in\mathfrak{S}_p} \mathbb{P}\left(\Gamma=c|Z^{(1)},\ldots,Z^{(n)}\right)
```

one can reasonably assume the data $`Z`$ was drawn from the model
``` math
\mathrm{N}_p(0,\pi_{c_{\text{map}}}(S))
```

where $`S = \frac{1}{n} \sum_{i=1}^n Z^{(i)}\cdot {Z^{(i)}}^\top`$

In such a case, we call $`c_{\text{map}}`$ the Maximum A Posteri (MAP).

#### Finding the MAP Estimator

The space of all permutations is enormous for bigger $`p`$ ($`p\ge9`$).
It is more reasonable to estimate the MAP in such a big space than to
calculate it precisely.

Metropolis-Hastings algorithm suggested by the authors of
[\[1\]](https://arxiv.org/abs/2004.03503) is a natural way to do it. To
see the discussion on it and other options available in `gips`, see
[`vignette("Optimizers")`](https://przechoj.github.io/gips/articles/Optimizers.md)
or its [pkgdown
page](https://przechoj.github.io/gips/articles/Optimizers.html).

#### Example

``` r

# Prepare model, multivariate normal distribution
p <- 6
number_of_observations <- 4
mu <- numeric(p)
sigma_matrix <- matrix(
  data = c(
    1.05, 0.8, 0.6, 0.4, 0.6, 0.8,
    0.8, 1.05, 0.8, 0.6, 0.4, 0.6,
    0.6, 0.8, 1.05, 0.8, 0.6, 0.4,
    0.4, 0.6, 0.8, 1.05, 0.8, 0.6,
    0.6, 0.4, 0.6, 0.8, 1.05, 0.8,
    0.8, 0.6, 0.4, 0.6, 0.8, 1.05
  ),
  nrow = p, byrow = TRUE
) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5,6)

# Generate example data from a model:
withr::with_seed(1234,{
  Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
})
# End of prepare model
```

Show/hide data preparation

Let’s say we have this data, `Z`. It has only $`4`$ observations, and
its dimension is $`p=6`$. Let’s assume `Z` was drawn from the normal
distribution with the mean $`(0,0,0,0,0,0)`$. We want to estimate the
covariance matrix:

``` r

dim(Z)
#> [1] 4 6
number_of_observations <- nrow(Z) # 4
p <- ncol(Z) # 6

# Calculate the covariance matrix from the data (assume the mean is 0):
S <- (t(Z) %*% Z) / number_of_observations

# Make the gips object out of data:
g <- gips(S, number_of_observations, was_mean_estimated = FALSE)

g_map <- find_MAP(g, optimizer = "brute_force")
#> ================================================================================
print(g_map)
#> The permutation (1,2,3,4,5,6)
#>  - was found after 720 log_posteriori calculations
#>  - is 157.903480698964 times more likely than the starting, () permutation.

S_projected <- project_matrix(S, g_map[[1]])
```

![](Theory_files/figure-html/example2_readme3-1.png)

We see the posterior probability
[\[1,(30)\]](https://arxiv.org/abs/2004.03503) has the biggest value for
the permutation $`(1,2,3,4,5,6)`$. It was over 150 times bigger than for
the trivial $`\text{id} = (1)(2)\ldots(p)`$ permutation. We interpret
that under the assumptions (centered Gaussian), it is over 150 times
more reasonable to assume the data `Z` was drawn from model
$`\mathrm{N}_p(0,\text{S_projected})`$ than from model
$`\mathrm{N}_p(0,\text{S})`$.

## References

\[1\] Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam.
“Model selection in the space of Gaussian models invariant by symmetry.”
The Annals of Statistics, 50(3) 1747-1774 June 2022. [arXiv
link](https://arxiv.org/abs/2004.03503); [DOI:
10.1214/22-AOS2174](https://doi.org/10.1214/22-AOS2174)
