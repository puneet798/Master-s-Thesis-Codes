# Master's-Thesis-Codes
This includes all the code used to plot the uncertainty in estimator of relative errors for a wide range of Ax=b problems.

Main script is plotLBUR to plot the figures. Sequence of the file run is plotLBUR -> UncerRatioPlotter -> crazyA, gmres11 and bicgtest.

crazyA generates a matrix (mostly non-symmetric positive definite) with a defined condition number and size of the matrix.

gmres11 and bicgtest are the codes for running GMRES and BiCG algorithm respectively.
