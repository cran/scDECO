# scDECO 0.1.1

* scdeco.cop: Added optional `offset1`, `offset2` arguments, which are added to the linear predictors for `mu1` and `mu2` (on the link scale).
* scdeco.cop: Changed starting values for `beta` parameters from all 0 to values informed by `mean(y)`.
* scdeco.cop: Changed starting values for `eta` parameters from all 0 to values informed by `mean(y == 0)`.
* scdeco.cop: Fixed issue where the `burn` and `thin` parameters were not actually used.
