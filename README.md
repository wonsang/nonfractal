# nonfractal: MATLAB Toolbox for estimating nonfractal connectivity 

## Description 
It is a MATLAB toolbox for estimating nonfractal connectivity from a set of time series such as resting state fMRI BOLD signals.

## Depends 
* Statistics Toolbox
* wmtsa

## Author 
Wonsang You (wsgyou@gmail.com)

## Getting Started
This toolbox was designed to work from MATLAB version 7.7 (R2008b) to the latest version with the Statistics toolbox. You can install this toolbox as follows.
1. Place all files of this toolbox at any directory, and add all sub-directories to the MATLAB paths.
2. The wmtsa toolbox should be installed a priori. It can be downloaded from http://www.atmos.washington.edu/~wmtsa/ for free.
3. Try to run the demo script "bfn_demo_nonfractal". This demo shows an example of estimating nonfractal connectivity of resting state fMRI BOLD signals of the rat brain using the maximum likelihood estimator.

## Functions
Let X be a NxQ matrix of Q time series with length N. It might be a set of BOLD signals corresponding to multiple ROIs of the brain which was extracted from DICOM images. Then, the function "bfn_mfin_ml" estimates the Hurst exponent and covariance matrix of short memory in multivariate time series using the maximum likelihood method as follows.

[H, nfcon] = bfn_mfin_ml(X);

H is the Hurst exponent, and nfcon is the nonfractal connectivity which is identical to the correlation matrix of short memory. On the other hand, the function "bfn_mfin_lms" estimates the nonfractal connectivity using the least-mean squares (LMS) method.

[H, nfcon] = bfn_mfin_lms(X);

For more options, run the commands "help bfn_mfin_ml" or "help bfn_mfin_lms".


## Reference
Please cite the following paper when using this toolbox.

Wonsang You, Sophie Achard, Joerg Stadler, Bernd Bruekner, and Udo
Seiffert, "Fractal analysis of resting state functional connectivity of
the brain," in 2012 International Joint Conference on Neural Networks,
2012. 

https://sites.google.com/site/wsgyou/publications/ijcnn2012

## Information
Please visit the webpage for more information:
https://sites.google.com/site/wsgyou/research/restingstate
