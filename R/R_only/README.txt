///////////////////////////////////////////////////////////////////////////////////
//Cascade Size Distributions: Why They Matter and How to Compute Them Efficiently//
/////////////////////////////////////////////////////////////////////////////////

The following code accompanies the paper with the above title.
In this folder, we provide an implementation of the main functions to compute cascade size distributions and conditional activation probabilities by SDP and TDA for the Independent Cascade Model (ICM).
It is written purely in R for readability and easy extensions.
Please note that we also provide more efficient C++ code, which can be called from R. We highly recommend C++ for large networks.

1) The file main.R contains a short tutorial how to use the main functions. 
2) The script randomGraphs.R provides functions for random graph generation, specifically, random trees, configuration model random graphs with scale free degree distribution and specified average degree, and random graphs with high clustering coefficient.
3) The file cascadeSizeDistr.R includes the main algorithms: SDP, BP, and TDA.

 