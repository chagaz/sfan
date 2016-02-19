# sfan
Selecting features as (network) nodes.

# Disclaimer
This work is still ongoing, and some parts have not been experimentally validated nor (for obvious reasons) peer-reviewed yet.

# Overview
The goal of this project is to provide a tool for various flavors of network-guided feature selection with graph cuts. This tool is meant for feature selection on data with the following characteristics:
- a relevance score, which characterizes the importance of each feature, is available. The underlying premise is that the relevance of a set of features is the sum of the relevances of that set. Some examples of relevance scores include Pearson's correlation and SKAT (Wu et al., 2011);
- a network (graph) is given over the features. In this network, nodes (vertices) represent features, while links (edges) connect two features that "work together" or "should tend to be selected together". The network can be weighted.

In this context, feature selection is done by finding the set of features that maximizing an objective function that is the sum of 
- the relevance of the selected features;
- a regularization term that (1) encourages sparsity (few features are selected) and (2) encourages the selected features to be connected to each other in the network.

This approach was first described in Azencott et al. (2013).

There are also several multi-task variants of this, in which one attempts to solve the same problem on multiple related tasks simultaneously. The features aren't necessarily the same for all tasks (although a large overlap is required for the approach to make sense) and the networks can also be task specific.

The first of these multitask approaches was described in Sugiyama et al. (2014).

These approaches have been developed with biomarker discovery in mind. In this case, the features are SNPs or other biomarkers (genes, epigenetic markers...) and we make use of biological networks. 

# Contributors
This work was initiated at the Centre for Computational Biology (CBIO) of MINES ParisTech and Institut Curie by Chloé-Agathe Azencott. Jean-Daniel Granet and Killian Poulaud worked on previous, private versions of this code. 
Please contact Chloé-Agathe Azencott at chloe-agathe.azencott@mines-paristech.fr for all inquiries.

# Requirements

# Installation

# Usage

# References
Azencott, C.-A., Grimm, D., Sugiyama, M., Kawahara, Y., and Borgwardt, K.M. (2013). Efficient network-guided multi-locus association mapping with graph cuts. Bioinformatics 29, i171–i179.

Sugiyama, M., Azencott, C., Grimm, D., Kawahara, Y., and Borgwardt, K. (2014). Multi-Task Feature Selection on Multiple Networks via Maximum Flows. In Proceedings of the 2014 SIAM International Conference on Data Mining, pp. 199–207.

Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011). Rare-Variant Association Testing for Sequencing Data with the Sequence Kernel Association Test. The American Journal of Human Genetics 89, 82–93.


