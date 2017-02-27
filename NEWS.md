# addreg 3.0 _(in progress)_
* `method` option added: fit using combinatorial EM or single parameter expanded EM algorithm
* `accelerate` option added: use `turboEM` to allow acceleration of EM algorithm
* Corrected an error in the calculation of deviance residuals for `family = negbin1`
* Corrected an error in the calculation of null deviance for `family = negbin1`
* (Minor change) Change `iter` to return a vector of two values

# addreg 2.0
* Added additive negative binomial regression (`family = negbin1`)
* (Minor change) Removed `bin.identity` link function; avoided CRAN error in a simpler way
* (Minor change) Documentation changes, including updated references

# addreg 1.2
* Introduced `bin.identity` link function: identical to the standard `identity` link, but does not cause an error when used with the `binomial` family.
* (Minor change) Corrections to avoid `codetools` errors.

# addreg 1.1
* The object returned by `nnpois()` now includes component `standard`.
* (Minor change) Various very minor corrections.

# addreg 1.0
* Added `na.action` argument to `addreg()` and `addreg.smooth()`