# gips

gips - Gaussian model Invariant by Permutation Symmetry

`gips` is an R package that allows finding the permutation symmetry group such that the variance of the data is invariant under it. Knowledge of such a permutation can drastically decrease the number of parameters needed to fit the model. That means it is possible to find the Gaussian model with more parameters than the number of observations with `gips`.


## `gips` will help you with two things:
1. Exploratory Data Analysis - with `gips`, you can find the permutation of features that does not change the covariance matrix.
2. Modeling - with `gips`, you can use the found permutation to fit the model accurately.

# Creddits

It will be developped by Przemysław Chojecki and Paweł Morgen under the leadership of Ph.D. Bartosz Kołodziejek.
