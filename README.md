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

The resulting optimization problem can be rewritten as a minimum cut problem. We solve it using a version of the Shekhostov implementation (Shekhovtsov & Hlavac, 2011) of the  Goldberg-Tarjan-Cherkassky algorithm (Cherkassky et al., 1995.)

This approach was first described in Azencott et al. (2013).

There are also several multi-task variants of this, in which one attempts to solve the same problem on multiple related tasks simultaneously. The features aren't necessarily the same for all tasks (although a large overlap is required for the approach to make sense) and the networks can also be task specific.

The first of these multitask approaches was described in Sugiyama et al. (2014).

These approaches have been developed with biomarker discovery in mind. In this case, the features are SNPs or other biomarkers (genes, epigenetic markers...) and we make use of biological networks. 

# Contributors
This work was initiated at the Centre for Computational Biology (CBIO) of MINES ParisTech and Institut Curie by Chloé-Agathe Azencott. Jean-Daniel Granet and Killian Poulaud worked on previous, private versions of this code. 
Please contact Chloé-Agathe Azencott at chloe-agathe.azencott@mines-paristech.fr for all inquiries.

# Requirements
* [cython](http://cython.org/)
* g++
* Python2.7
* Python header files (libpython-dev)
* Python libraries/packages/ecosystems:
  * [NumPy](http://www.numpy.org/)
  * [SciPy](http://www.scipy.org/)
  * [PyTables](http://www.pytables.org/), at least for data generation.

# Installation
```
cd code
./__build_gt_maxflow.sh
```
Create `gt_maxflow.so`.

# Testing
For testing (doctest) the core optimization module run by `code/multitask_sfan.py`:
```
cd code/multitask_sfan.py
python multitask_sfan.py -t
```

# Usage
## Core optimization
The core optimization (for given regularization parameters) is run by `code/multitask_sfan.py`. See `code/test_multitask_sfan.sh` for usage.
Example:
```
 python multitask_sfan.py --num_tasks 2 --networks ../data/simu_01/simu_01.network.dimacs \
       --node_weights ../data/simu_01/simu_01.scores_0.txt ../data/simu_01/simu_01.scores_1.txt \
       --correlation_matrix ../data/simu_01/simu_01.task_similarities.txt -l 0.001 -e 0.02 -m 0.01
```

## Data generation
`generate_data.py` generates synthetic data for experiments:
* a modular network (modules are fully connected) over the features;
* a genotype (SNP) matrix X of random integers between 0 and 2;
* a similarity matrix Omega between tasks;
* causal features (SNPs), with corresponding weights, generated so as to respect the covariance structure given by
Omega;
* the corresponding k phenotypes and vectors of node weights (computed as Pearson's correlation).

Example:
```
 python generate_data.py -k 2 -m 1000 -n 10 ../data/simu_02 simu_02
```

# File formats
## Relevance scores
Each line is the (floating point) relevance score of the corresponding node (in the same order as for the network file).

Example: `data/simu_multitask_01.scores_0.txt`.

## Networks
Networks built over features are given as "DIMACS for maxmimum flow" files. 
For details about this format see [DIMACS maximum flow problems](http://lpsolve.sourceforge.net/5.5/DIMACS_maxf.htm).

Network files begin with a header of the form:
```
p max number_of_nodes number_of_edges
```
and each arc is listed as
```
a node_a node_b edge_weight
```
__Important__
* An undirected edge must be listed twice, first as `a node_a node_b edge_weight`, then as `a node_b node_a edge_weight`.
* Because of the maxflow implementation, edges must be __ordered__ by increasing values of `node_a` and then of `node_b`.
* Node indices start at 1.

Example: `data/simu_multitask_01.network.dimacs`. Note that edge weights can be given as floats.

We use this format because it is the one that `gt_maxflow` (Shekhovtsov & Hlavac, 2011) understands. 

## Correlation between tasks
number_of_tasks x number_of_tasks symmetric matrix (space-separated columns).

Example: `data/simu_multitask_01.task_similarities.txt`. 

## Output of `code/multitask_sfan.py`
Commented lines give the parameter values.
Then each line corresponds to one tasks, and is a space-separated list of node indices, __starting at 1__.

# References
Azencott, C.-A., Grimm, D., Sugiyama, M., Kawahara, Y., and Borgwardt, K.M. (2013). Efficient network-guided multi-locus association mapping with graph cuts. Bioinformatics 29, i171–i179.

Cherkassky, Boris V., and Andrew V. Goldberg. “On Implementing Push-Relabel Method for the Maximum Flow Problem.” In Integer Programming and Combinatorial Optimization, edited by Egon Balas and Jens Clausen, 920:157–71. Berlin, Heidelberg: Springer Berlin Heidelberg, 1995. 

Shekhovtsov A. and Hlavac V. (2011). A Distributed Mincut/Maxflow Algorithm Combining Path Augmentation and Push-Relabel. EMMCVPR 2011 / Technical Report CTU--CMP--2011--03. http://cmp.felk.cvut.cz/~shekhovt/d_maxflow

Sugiyama, M., Azencott, C., Grimm, D., Kawahara, Y., and Borgwardt, K. (2014). Multi-Task Feature Selection on Multiple Networks via Maximum Flows. In Proceedings of the 2014 SIAM International Conference on Data Mining, pp. 199–207.

Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011). Rare-Variant Association Testing for Sequencing Data with the Sequence Kernel Association Test. The American Journal of Human Genetics 89, 82–93.


