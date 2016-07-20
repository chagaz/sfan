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

Technical details are available under `paper/tech_notes.pdf` (to be compiled from `paper/tech_notes.tex`) 

# Contributors
This work was initiated at the Centre for Computational Biology (CBIO) of MINES ParisTech and Institut Curie by Chloé-Agathe Azencott. Jean-Daniel Granet and Killian Poulaud worked on previous, private versions of this code. 
Please contact Chloé-Agathe Azencott at chloe-agathe.azencott@mines-paristech.fr for all inquiries.

# Requirements
* [cython](http://cython.org/)
* g++ (code tested with gcc4.8.4)
* Python2.7
* Python header files (libpython-dev)
* Python libraries/packages/ecosystems:
  
  * [NumPy](http://www.numpy.org/)
  * [SciPy](http://www.scipy.org/)
  * [PyTables](http://www.pytables.org/), at least for data generation.
  * [scikit-learn](http://scikit-learn.org/stable)
  * [Matplotlib](http://matplotlib.org/)
    
  From standard : 
  * [argparse](https://docs.python.org/2/library/argparse.html)
  * [doctest](https://docs.python.org/2/library/doctest.html)
  * [glob](https://docs.python.org/2.7/library/glob.html)
  * [logging](https://docs.python.org/2.7/library/logging.html)
  * [math](https://docs.python.org/2.7/library/math.html)
  * [operator](https://docs.python.org/2.7/library/operator.html)
  * [os](https://docs.python.org/2.7/library/os.html)
  * [random](https://docs.python.org/2.7/library/random.html)
  * [shlex](https://docs.python.org/2.7/library/shlex.html)
  * [shutil](https://docs.python.org/2.7/library/shutil.html)
  * [subprocess](https://docs.python.org/2.7/library/subprocess.html)
  * [sys](https://docs.python.org/2.7/library/sys.html)
  * [tempfile](https://docs.python.org/2.7/library/tempfile.html)
  * [time](https://docs.python.org/2.7/library/time.html)


* python-dev
# Installation
```
cd code
./__build_gt_maxflow.sh
```
Create `gt_maxflow.so`.

Modifications to ```__build_gt_maxflow.sh```:
* You might need to replace   ```-I/usr/include/python2.7``` with the path to your ```Python.h``` (Python development header).
* On some architectures, or if you have an older version of gcc, you will need to add the flag ```-D__USE_XOPEN2K8``` to the ```g++``` line. This is linked to gcc's fixinclude mechanism, and usually manifests itself as errors regarding the type name ```__locale_t```.
* If you're getting ```undefined reference to 'clock_gettime' and 'clock_settime'``` errors, you might need to add the flag ```-lrt``` to the ```g++``` line.

# Testing
For testing (doctest) the core optimization module run by `code/multitask_sfan.py`:
```
cd code
python multitask_sfan.py -t -v
```

# Usage
## Core optimization
The core optimization (for given regularization parameters) is run by `code/multitask_sfan.py`. See `code/test_multitask_sfan.sh` for usage.

Single task example:
```
 python multitask_sfan.py --num_tasks 1 --networks ../data/simu_01/simu_01.network.dimacs \
       --node_weights ../data/simu_01/simu_01.scores_0.txt -l 0.001 -e 0.02 --output runtime.txt
```

Multitask example:
```
 python multitask_sfan.py --num_tasks 2 --networks ../data/simu_01/simu_01.network.dimacs \
       --node_weights ../data/simu_01/simu_01.scores_0.txt ../data/simu_01/simu_01.scores_1.txt \
       --covariance_matrix ../data/simu_01/simu_01.task_similarities.txt -l 0.001 -e 0.02 -m 0.01
```
The user can either provide a covariance or a precision matrix between tasks. The covariance matrix encodes a notion of similarity between the tasks. The precision matrix is its inverse, and its off-diagonal entries can be interpreted as the normalized opposite of the partial correlation between the corresponding tasks. The methods only makes it possible to account for positive or non-existant partial correlations, meaning that positive off-diagonal entries of the precision matrix, if any, will be thresholded to 0.

If no covariance nor precision matrix is given, a precision matrix with (<number of tasks>-1+epsilon) on the diagonal and -1 off the diagonal is used and the value of eta is adjusted to match the formulation of MultiSConES by Sugiyama et al. (2014).

## Data generation
`code/generate_data.py` generates synthetic data for experiments:
* a modular network (modules are fully connected) over the features;
* a genotype (SNP) matrix X of random integers between 0 and 2;
* a covariance matrix Omega of similarities between tasks;
* causal features (SNPs), with corresponding weights, generated so as to respect the covariance structure given by
Omega (inverse of the precision matrix);
* the corresponding k phenotypes and vectors of node weights (computed as Pearson's correlation).

Example:
```
cd code
python generate_data.py -k 3 -m 1000 -n 50 ../data/simu_synth_01 simu_01 --verbose
```

## Synthetic data experiments -- TO BE COMPLETED
`code/evaluation_framework.py` contains methods and classes needed for evaluation (determining cross-validation sets, computing performance, etc.)

`code/synthetic_data_experiment.py` runs experiments on synthetic data.

`code/handle-output.py` uses output of `code/synthetic_data_experiment.py` to produce results 
- For each fold : 
    - Fit a ridge regression model between finally selected features genotype of train set and their phenotype
    - Predict phenotype of test set using  this model with finally selected features genotype of test set
- Output values of each criterion for each fold of each plot 
- Output summary 
- Boxplot results

Example : 
```
cd code
python synthetic_data_experiments.py -k 3 -m 200 -n 100 -r 10 -f 10 -s 10 \
             ../data/simu_synth_01 ../results/simu_synth_01 simu_01 --verbose
python handle-output.py -k 3 -m 200 -n 100 -r 10 -f 10 -s 10 \
             ../data/simu_synth_01 ../results/simu_synth_01 simu_01 --verbose
```

### Usage on SGE cluster 
Some nodes of the SGE cluster of CBIO has problems using PyTables, so generate data before running experimentation and ensure DATA_GEN flag of ```synthetic_data_experiment.py``` is ```False```. Moreover, ensure ```SEQ_MODE``` flag is ```False``` and set a ```tmp_dir``` .

```synthetic_data_experiment.py``` will use a qsub job for each fold. 

#### Runtime performances (on SGE cluster only) 
Ensure TIME_EXP flag at the begining of ```synthetic_data_experiment.py``` is True. This qsub job on node 15-24 which have the same spec.



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

## Covariance (similarity) between tasks
number_of_tasks x number_of_tasks symmetric matrix (space-separated columns).

Example: `data/simu_multitask_01.task_similarities.txt`. 

## Output of `code/multitask_sfan.py`
Commented lines give the parameter values.
Then each line corresponds to one tasks, and is a space-separated list of node indices, __starting at 1__.

## Output of 'code/synthetic_data_experiments.py'


### For each algo ('sfan', 'msfan_np', 'msfan'):
* <simu_id>.<algo>.maxRSS : 
  maxRSS used during feature selection using all training set and optimal parameters.
  One value per line, one line per fold, all repeats mixed.

* <simu_id>.<algo>.timing : 
  Timing infos when runing feature selection using all training set and optimal paramters :
  ```
  Task (<task_idx>) computation time : <value>
  Task average computation time: <value>
  Standard deviation computation time: <value>
  Network building time: <value>
  gt_maxflow computation time: <value>
  ```
  for each fold. 



#### For each repeat : 

##### For each fold : 
* <repeat_idx>/<simu_id>.<algo>.fold_<fold_idx>.parameters : 
list of optimal parameters retain for the fold
* <repeat_idx>/<simu_id>.<algo>.fold_<fold_idx>.selected_features : 
space separated list of features
one line per task 
* <repeat_idx>/<simu_id>.<algo>.fold_<fold_idx>.ss.maxRSS : 
One value of max RSS per line, one line per subsample
* <repeat_idx>/<simu_id>.<algo>.fold_<fold_idx>.ss.process_time : 
One value of max RSS per line, one line per subsample


###### For each measure (pvv, tpr) :
<repeat_idx>/<simu_id>.<algo>.fold_<fold_idx>.<measure> : 
space-separated list of value measure
* <repeat_idx>/<simu_id>.<algo>.fold_<fold_idx>.ppv
* <repeat_idx>/<simu_id>.<algo>.fold_<fold_idx>.tpr
  




## Output of 'code/handle-output.py'

### Results are saved in following files : 
* <simu_id>.results :
  Average/standard deviation values for: consistency index, RMSE, PPV and TPR, as LaTeX table.



#### For each algo ('sfan', 'msfan_np', 'msfan'):

* <simu_id>.<algo>.rmse : 
  List of final RMSEs 
  one line per repeat
  for each repeat, one value per task
* <simu_id>.<algo>.consistency : 
  List of final Consistency Indices 
  one line per repeat
  for each repeat, one value per task

#### For each classification measure (accuracy (acc), Mathieu coefficient (mcc), Prositive Predictive Value (ppv) and True Positive Value (tpr) )
a file named <simu_id>.<algo>.<measure> :
Space-separated lists of PPVs (one value per task and per fold),
each line corresponds to one repeat. 
* <simu_id>.<algo>.acc
* <simu_id>.<algo>.mcc
* <simu_id>.<algo>.ppv
* <simu_id>.<algo>.tpr 


#### For each repeat, fold and task : 
* <repeat_idx>/<simu_id>.<algo>.fold_<fold_idx>.task_<task_idx>.predicted : 
phenotype prediction of the test set using a ridge-regression trained with the selected features only.
One value per line, one line per sample

    
    
    
### Charts :
For each measure (ci, mcc, ppv, tpr, rmse, acc), a chart named <simu_id>.<measure>.png : one boxplot per algorithm, grouped per task, with error bars. 
* <simu_id>.ci.png
* <simu_id>.rmse.png
* <simu_id>.acc.png
* <simu_id>.mcc.png
* <simu_id>.ppv.png
* <simu_id>.tpr.png


# References
Azencott, C.-A., Grimm, D., Sugiyama, M., Kawahara, Y., and Borgwardt, K.M. (2013). Efficient network-guided multi-locus association mapping with graph cuts. Bioinformatics 29, i171–i179.

Cherkassky, Boris V., and Andrew V. Goldberg. “On Implementing Push-Relabel Method for the Maximum Flow Problem.” In Integer Programming and Combinatorial Optimization, edited by Egon Balas and Jens Clausen, 920:157–71. Berlin, Heidelberg: Springer Berlin Heidelberg, 1995. 

Shekhovtsov A. and Hlavac V. (2011). A Distributed Mincut/Maxflow Algorithm Combining Path Augmentation and Push-Relabel. EMMCVPR 2011 / Technical Report CTU--CMP--2011--03. http://cmp.felk.cvut.cz/~shekhovt/d_maxflow

Sugiyama, M., Azencott, C., Grimm, D., Kawahara, Y., and Borgwardt, K. (2014). Multi-Task Feature Selection on Multiple Networks via Maximum Flows. In Proceedings of the 2014 SIAM International Conference on Data Mining, pp. 199–207.

Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011). Rare-Variant Association Testing for Sequencing Data with the Sequence Kernel Association Test. The American Journal of Human Genetics 89, 82–93.


