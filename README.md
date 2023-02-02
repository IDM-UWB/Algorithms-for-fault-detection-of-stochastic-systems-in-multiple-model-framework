# Algorithms for fault detection of stochastic systems in multiple model framework
Variable Structure (VS) Interacting Multiple Model (IMM) state estimator implementation

Based on series of papers by X. Rong Li. The main source is __X. Rong Li (200) Chapter 10 Engineer's Guide to Variable-Structure Multiple-Model Estimation for Tracking__.

## The implementation of VSIMM state estimator is illustrated in two examples
* _example_maneuveringtarget.m_ illustrates the use of the VS IMM for state estimation of a maneuvering target that undergoes different modes of motions (straight, turning with different turn rates). Details can be found in X. Rong Li (1999) Multiple-Model Estimation with Variable Structure Part iV - Design and Evaluation of Model-Group Switching Algorithm, IEEE Transactions on Aerospace and Electronic Systems, vol. 35, no. 1, pp. 242-254.
* _example_faultdetection.m_ illustrates the use of the VS IMM for fault detection of additive faults described in multiple model framework. Details can be found in I. Puncochar and O. Straka (2022) Parity-Space and Multiple-Model based Approaches to Measurement Fault Estimation, In Proceeding of the 16th European Workshop on Advanced Control and Diagnosis, Nancy, France.

## Main component of the VSIMM
* _vsimmkff.m_  implements the filtering step of the VS IMM using the Kalman filter
* _vsimmkfp.m_ implements the prediction step of the VS IMM using the Kalman filter
* _vsimmadatpre.m_ implements the adaptation of set of active models (extension)
* _vsimmadaptpost.m_ implements the adaptation of set of active models (discarding)
* _kfp.m_ implements the predictive step of the Kalman filter
* _kff.m_ implements the filtering step of the Kalman filter
* _normalizeweightexp.m_ implements numerically robust normalization of weights to probabilities


## Auxiliary functions
* _gendrnd.m_ implements the random number generator from discrete distribution with given probabilities
* _normrndm.m_ implements the vectorized version of standard normrnd.m
* _isnotempty.m_ implements the test for a nonempty array
* _tilefigure.m_ implements the tiling figures
* _myprint.m_ implements the printing figures into files
* _createmultiplemodels.m_ implements creation of models for the second example
* _generateparametervariants.m_ implements the construction of multiple model structure for given discrete set values of parameters of the model
* _mmmerge.m_ implements merging of mean values, covariance matrices and probabilities for terminal subsequences of the multiple hypothesis tree
* _ps2pfs.m_ implements merging of probabilities for terminal subsequences of the multiple hypothesis tree
* _smmkff.m_ implements the filtering step of the Kalman filter for switch multiple model

