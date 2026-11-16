# XYOnion: Layer-Based Splitting into Calibration and Validation Subsets

This repository contains the MATLAB code developed for the study *XYOnion: a layer-based method for splitting datasets into calibration and validation subsets*. The method combines the onion-layer approach with distances that jointly consider predictor space (`X`) and response space (`y`).

## Method overview

XYOnion assigns samples to calibration and validation subsets through successive layers of the data space. This implementation uses the DISTSLCT algorithm to reduce the computational cost of sample selection, avoiding calculation of all inter-sample distances at each iteration. Remaining samples are assigned randomly according to the requested calibration fraction.

The repository includes:

- `XYonion.m`: MATLAB function, including its sample-selection routines.
- The associated research article in PDF format.

## Purpose

The method supports the construction of representative calibration and validation subsets for chemometric modelling by considering both predictors and reference values during dataset splitting.

## Usage

Add the repository folder to the MATLAB path and run:

```matlab
rng(1); % For reproducibility of the random assignment
split = XYonion(X, y, 0.66);

Xcal = X(split, :);
ycal = y(split, :);
Xval = X(~split, :);
yval = y(~split, :);
```

The full function signature is:

```matlab
split = XYonion(X, y, fraction, nonion, mahalanobis);
```

- `X`: numeric matrix with samples in rows and variables in columns.
- `y`: column vector of reference values corresponding to the rows of `X`.
- `fraction`: target fraction of samples assigned to calibration; default: `0.66`.
- `nonion`: number of onion layers; default: `3`.
- `mahalanobis`: set to `1` to whiten `X` using its covariance matrix before sample selection; default: `0`.

The output `split` is a logical column vector: `true` indicates calibration samples and `false` indicates validation samples. The final subset sizes may differ slightly from the requested proportion due to rounding. The optional covariance whitening requires a nonsingular covariance matrix.

## Citation

If you use this code, please cite:

> Ezenarro, J. (2025). *XYOnion: a layer-based method for splitting datasets into calibration and validation subsets*. Analytica Chimica Acta, 1364, 344229. https://doi.org/10.1016/j.aca.2025.344229

## License

The MATLAB code is dedicated to the public domain under [CC0 1.0 Universal](LICENSE.txt). Reuse, modification, and redistribution are permitted for any purpose, including commercially. Citation of the associated article is appreciated but is not a requirement of CC0.

The included research article is published under CC BY-NC-ND 4.0, as stated in the PDF, and is not covered by the code's CC0 dedication.
