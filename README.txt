=============================================================================================
MATLAB code for: 
XYOnion — a layer-based method for splitting datasets into calibration and validation subsets
=============================================================================================

XYOnion splits a dataset into calibration and test sets using an onion-layer
idea guided by SPXY distances (joint consideration of predictor space `X` and
response `y`).

This version uses the DISTSLCT algorithm to avoid calculating all inter-sample
distances in each iteration.

Function signature
------------------

    split = XYonion(X, y, fraction, nonion, mahalanobis)

Inputs
------
- `X` — `m × p` double. Samples in rows, variables in columns.
- `y` — `m × 1` (or `1 × m`) double. Reference values aligned with rows of `X`.
- `fraction` — scalar in (0,1). Target fraction of samples assigned to calibration.
- `nonion` *(optional)* — positive integer. Number of onion layers. Default: 3.
- `mahalanobis` *(optional)* — logical/scalar {0,1}. If 1, whiten X so Euclidean
  distances act like Mahalanobis. Default: 0.

Output
------
- `split` — `m × 1` logical. `true` (or 1) = calibration, `false` (or 0) = test.


Reference and citation
----------------------
Author: Jokin Ezenarro  
Version: v1.0 (10/6/2025)  
Ref.: https://doi.org/10.1016/j.aca.2025.344229
