# Algorithms for fault detection of stochastic systems in multiple model framework
Variable Structure (VS) Interacting Multiple Model (IMM) state estimator implementation

Based on series of papers by X. Rong Li. The main source is _X. Rong Li (200) Chapter 10 Engineer's Guide to Variable-Structure Multiple-Model Estimation for Tracking_.

## The implementation of VSIMM state estimator is illustrated in two examples
* _example_maneuveringtarget.m_ - illustrates the use of the VS IMM for state estimation of a maneuvering target that undergoes different modes of motions (straight, turning with different turn rates). Details can be found in X. Rong Li (1999) Multiple-Model Estimation with Variable Structure Part iV - Design and Evaluation of Model-Group Switching Algorithm, IEEE Transactions on Aerospace and Electronic Systems, vol. 35, no. 1, pp. 242-254.
* _example_faultdetection.m_ - illustrates the use of the VS IMM for fault detection of additive faults described in multiple model framework. Details can be found in I. Puncochar and O. Straka (2022) Parity-Space and Multiple-Model based Approaches to Measurement Fault Estimation, In Proceeding of the 16th European Workshop on Advanced Control and Diagnosis, Nancy, France.

## Main component of the VSIMM
* vsimmkff.m - implements the filtering step of the VS IMM using the Kalman filter
* vsimmkfp.m - implements the prediction step of the VS IMM using the Kalman filter
* vsimmadatpre.m - implements the adaptation of set of active models (extension)
* vsimmadaptpost.m - implements the adaptation of set of active models (discarding)
* kfp.m - implements the predictive step of the Kalman filter
* kff.m - implements the filtering step of the Kalman filter
* normalizeweightexp.m - implements numerically robust normalization of weights to probabilities


## Auxiliary functions
* gendrnd.m - implements the random number generator from discrete distribution with given probabilities
* normrndm.m - implements the vectorized version of standard normrnd.m
* isnotempty.m - implements the test for a nonempty array
* tilefigure - implements the tiling figures
* myprint - implements the printing figures into files
* createmultiplemodels.m - implements creation of models for the second example
* generateparametervariants.m - implements the construction of multiple model structure for given discrete set values of parameters of the model
* mmmerge.m - implements merging of mean values, covariance matrices and probabilities for terminal subsequences of the multiple hypothesis tree
* ps2pfs.m - implements merging of probabilities for terminal subsequences of the multiple hypothesis tree
* smmkff - implements the filtering step of the Kalman filter for switch multiple model
