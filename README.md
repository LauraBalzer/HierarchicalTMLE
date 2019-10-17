# HierarchicalTMLE
R code to generate simulated data and implement the hierarchical TMLEs described in Balzer et al: "A New Approach to Hierarchical Data Analysis: Targeted Maximum Likelihood Estimation for the causal effect of a cluster-level exposure". 
Original posting June 8, 2017: https://arxiv.org/abs/1706.02675v1. Updates based on manuscript revision March 12, 2018. Revised manuscript available upon request. 

Published paper available at L.B. Balzer, W. Zheng, M.J. van der Laan, M.L. Petersen, and the SEARCH Collaboration. A new
approach to hierarchical data analysis: Targeted maximum likelihood estimation for the causal effect of a
cluster-level exposure. Statistical Methods in Medical Research, 28(6):1761–1780, 2018. https://journals.sagepub.com/doi/abs/10.1177/0962280218774936#

Warning: 
Be sure to check that the lines using the aggregate function on lines 179, 212, 585 are returning data.frames with the same number of columns as the original (non-aggregated) dataset. Depending on how your data are organized, this might be a particular issue for line 585.
