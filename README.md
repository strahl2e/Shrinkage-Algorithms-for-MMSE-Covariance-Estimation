As part of the coursework for Advanced Topics in Statistical Learning, Khaoula El Mekkaoui (https://github.com/Mekkhaoula) and I wrote this MATLAB code to replicate the experiments in a research paper on covariance (shinkage) estimation and to test the methods on real stock-market data - predictions were pretty good, but we're unfortunately not yet millionares, so wasn't quite that good :)

Prof. Esa Ollila used this work as a demo in future lectures, so it was quite good. :)

Experiments based on the paper Chen Y., Wiesel A., Eldar Y., Hero A.. Shrinkage Algorithms for MMSE Covariance Estimation. IEEE Transactions on Signal Processing, 2010.

The MATLAB files above include several stand-alone scripts for the experiments:

MSE_simulations.m - The mean-square error (MSE) simulations for covariance estimation from the paper.

PortfolioManShrinkCovEst.m - Applying methods to stock portfolio management, where investments can be balanced based on negative correlation estimates to reduce risk.


Application to financial portfolio managament as an example of the utility of these algorithms in a real-world setting.
