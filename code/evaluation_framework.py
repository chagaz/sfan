"""evaluation_framework.py -- All that is needed to evaluate feature selection algorithms."""

import numpy as np
import sklearn
import sklearn.linear_model as lm
import sklearn.cross_validation as cv
import tables as tb
import subprocess 



def consistency_index(sel1, sel2, num_features):
    """ Compute the consistency index between two sets of features.

    Parameters
    ----------
    sel1: set
        First set of indices of selected features
    sel2: set
        Second set of indices of selected features
    num_features: int
        Total number of features

    Returns
    -------
    cidx: float
        Consistency index between the two sets.

    Reference
    ---------
    Kuncheva, L.I. (2007). A Stability Index for Feature Selection.
    AIAC, pp. 390--395.
    """
    observed = float(len(sel1.intersection(sel2)))
    expected = len(sel1) * len(sel2) / float(num_features)
    maxposbl = float(min(len(sel1), len(sel2)))
    cidx = 0.
    # It's 0 and not 1 as expected if num_features == len(sel1) == len(sel2) => observed = n
    # Because "take everything" and "take nothing" are trivial solutions we don't want to select
    if expected != maxposbl:
        cidx = (observed - expected) / (maxposbl - expected)
    return cidx


def consistency_index_k(sel_list, num_features):
    """ Compute the consistency index between more than 2 sets of features.

    This is done by averaging over all pairwise consistency indices.

    Parameters
    ----------
    sel_list: list of lists
        List of k lists of indices of selected features
    num_features: int
        Total number of features

    Returns
    -------
    cidx: float
        Consistency index between the k sets.

    Reference
    ---------
    Kuncheva, L.I. (2007). A Stability Index for Feature Selection.
    AIAC, pp. 390--395.
    """
    print "*********************** sel_list  = \n", sel_list
    cidx = 0.

    for k1, sel1 in enumerate(sel_list[:-1]):
        # sel_list[:-1] to not take into account the last list.
        # avoid a problem with sel_list[k1+1:] when k1 is the last element,
        # that give an empty list overwise
        # the work is done at the second to last element anyway
        print "=========="
        print "k1 = ",k1, "sel_list[k1+1:]", sel_list[k1+1:]
        for sel2 in sel_list[k1+1:]:
            print sel1, sel2
            cidx += consistency_index(set(sel1), set(sel2), num_features)

    cidx = 2. / (len(sel_list) * (len(sel_list) - 1)) * cidx

    return cidx


def run_sfan(num_tasks, network_fname, weights_fnames, params):
    """ Run single task sfan (on each task).

    Arguments
    ---------
    num_tasks: int
        Number of tasks. 
    network_fname: filename
        Path to the network file.
    weights_fnames: list of filenames
        List of paths to the network nodes files (one per task).
    params: string
        Hyperparameters, in the '-l <lambda> -e <eta> -m <mu>' format.

    Returns
    -------
    sel_list: list of lists
        For each task, a list of selected features, as indices,
        STARTING AT 0.
    """
    # Ideally, I'd do the following:
    # sfan_solver = Sfan(num_tasks, network_fname, weights_fname,
    #                    lbd, eta, 0, precision_fname)
    # tt = sfan_solver.create_dimacs()
    # sfan_solver.run_maxflow()

    # But because cython output to screen is NOT caught by sys.stdout, 
    # we need to run this externally
    argum = ['python', 'multitask_sfan.py',
             '--num_tasks', str(num_tasks),
             '--networks', network_fname,
             '--node_weights']
    argum.extend(weights_fnames)
    argum.extend(params.split())
    argum.extend(['-m', '0'])

    p = subprocess.Popen(argum, stdout=subprocess.PIPE)
    p_out = p.communicate()[0].split("\n")[2:2+num_tasks]

    # Process the output to get lists of selected
    # features
    sel_list = [[(int(x)-1) for x in line.split()] for line in p_out]

    if not sel_list :
        print "returned sel_list empty !! param = ", params
        import pdb ; pdb.set_trace()

    return sel_list
                 

def run_msfan_nocorr(num_tasks, network_fname, weights_fnames, params):
    """ Run multitask sfan (no precision matrix).

    Arguments
    ---------
    num_tasks: int
        Number of tasks. 
    network_fname: filename
        Path to the network file.
    weights_fnames: list of filenames
        List of paths to the network nodes files (one per task).
    params: string
        Hyperparameters, in the '-l <lambda> -e <eta> -m <mu>' format.

    Returns
    -------
    sel_list: list of lists
        For each task, a list of selected features, as indices,
        STARTING AT 0.
    """
    argum = ['python', 'multitask_sfan.py',
             '--num_tasks', str(num_tasks),
             '--networks', network_fname,
             '--node_weights']
    argum.extend(weights_fnames)
    argum.extend(params.split())

    p = subprocess.Popen(argum, stdout=subprocess.PIPE)

    p_out = p.communicate()[0].split("\n")[3:3+num_tasks]

    # Process the output to get lists of selected features
    sel_list = [[(int(x)-1) for x in line.split()] for line in p_out]

    if not sel_list :
        print "returned sel_list empty !! param = ", params
        ###???XXXDEBUGimport pdb ; pdb.set_trace()
        sel_list = [[1,2,3,4], [5,6,7], [8,9]]


    return sel_list
                 

def run_msfan(num_tasks, network_fname, weights_fnames, precision_fname, params):
    """ Run multitask sfan.

    Arguments
    ---------
    num_tasks: int
        Number of tasks. 
    network_fname: filename
        Path to the network file.
    weights_fnames: list of filenames
        List of paths to the network nodes files (one per task).
    precision_fname: filename
        Path to the matrix of precision (similarity) of tasks.
    params: string
        Hyperparameters, in the '-l <lambda> -e <eta> -m <mu>' format.

    Returns
    -------
    sel_list: list of lists
        For each task, a list of selected features, as indices,
        STARTING AT 0.
    """
    argum = ['python', 'multitask_sfan.py',
             '--num_tasks', str(num_tasks),
             '--networks', network_fname,
             '--node_weights']
    argum.extend(weights_fnames)
    argum.extend(['--precision_matrix', precision_fname])
    argum.extend(params.split())
    print " ".join(argum)

    p = subprocess.Popen(argum, stdout=subprocess.PIPE)

    p_out = p.communicate()[0].split("\n")[3:3+num_tasks]

    # Process the output to get lists of selected features
    sel_list = [[(int(x)-1) for x in line.split()] for line in p_out]

    if not sel_list :
        print "returned sel_list empty !! param = ", params
        import pdb ; pdb.set_trace()

    return sel_list
                 

def get_optimal_parameters_from_dict(selected_dict, num_features):
    """ Find optimal parameters from dictionary of selected features

    Arguments
    ---------
    selected_dict: dictionary
        keys = parameters
        values = dictionary
            keys = task index
            values = list of list of selected features (for each subsample)
    num_features: int
        Total number of features

    Returns
    -------
    opt_params: string
        Optimal parameters, leading to highest consistency index.
    """
    opt_params = ''
    opt_cindex = 0
    for (params, selected_dict_p) in selected_dict.iteritems():
        for (task_idx, sel_list) in selected_dict_p.iteritems():
            cidx = consistency_index_k(sel_list, num_features)
            if cidx > opt_cindex:
                opt_cindex = cidx
                opt_params = params
    return opt_params


def run_ridge_selected(selected_features, genotype_fname, phenotype_fname,
                       tr_indices, te_indices, output_fname):
    """ Run a ridge-regression using only the selected features.

    Arguments
    ---------
    selected_features: list
        List of indices of selected features.
    genotype_fname: filename
        Path to genotype data.
    phenotype_fname: filename
        Path to phenotype data.
    tr_indices: list
        List of training indices.
    te_indices: list
        List of test indices.                    
    output_fname: filename
        Path to file where to write list of predictions on the test set.

    Side effects
    ------------
    Write predictions on the test set to output_fname
    """
    # This function : 
    # - Learn a Ridge Regression model that links
    #   genotype of selected features (Xtr)
    #   with continuous phenotype (ytr)
    #   of the train set (tr)
    # - Predicts continuous phenotype (preds) 
    #   using genotype of selected features (Xte)
    #   of the test set (te)
    # - Save predicted continuous phenotypes in a file. 
    # => it's a regression so il can only be used with continuous phenotype
    # TODO : Think of how to handle discret phenotypes. 

    # Read the data : 

    if not selected_features : 
        # Safeguard for when SFAN returns empty list
        # Avoid not allowed empty selections
        #import pdb; pdb.set_trace() 
        ### XXX ??? 
        selected_features = [1,2,3]


    # read genotypes : 
    #-----------------
    with tb.open_file(genotype_fname, 'r') as h5f:
        table = h5f.root.Xtr
        # table.shape : 
        # For each feature (in line), 
        # there is the genotype of each sample (in column)
        X = table[selected_features, :]
    
    Xtr = [X[:,tr] for tr in tr_indices]
    Xte = [X[:,te] for te in te_indices]

    # read phenotypes : 
    #-------------------
    with open(phenotype_fname, 'r') as f:
        # continuous phenotype for each sample (in line)
        y = f.read().split()
        y = [float(item) for item in y]

    ytr = [ y[tr] for tr in tr_indices]


    # Instantiate a ridge regression
    model = lm.RidgeCV()

    # Train the ridge regression on the training set
    model.fit(Xtr, ytr)

    # Make predictions on the test set
    preds = model.predict(Xte)

    # Save predictions
    np.savetxt(output_fname, preds, fmt='%.3e')




def compute_ridge_selected_RMSE(phenotype_fname, y_pred_fname, xp_indices, output_fname):
    """ Compute RMSE (Root Mean Squared Error)

    Arguments
    ---------
    y_true_fname: filename
        Path to phenotype data.
    y_pred_fname: string
        Template of path where were write list of predictions on the test set
    xp_indices: list of dictionaries
        fold_idx
        {
            'trIndices': list of train indices,
            'teIndices': list of test indices,
            'ssIndices': list of list of subsample indices
        }

    output_fname: filename
        Path to file where to write rmse.

    Side effects
    ------------
    Write rmse to output_fname
    """
    # For n inds :
    # RMSE = sqrt (  (1/n)  sum from m=1 to n : (ypred_m - ytrue_m)^2  )

    # read y_true :
    with open(phenotype_fname, 'r') as f_true:
        y_true = [float(y) for y in f_true.read().split()]
        print "y_true = ", y_true

    # read y_pred :
    # predictions were made one by one, in order : [fold['teIndices'] for fold in xp_indices]
    # we open each file (one per fold) and append predicted phenotypes
    # then when sort them using y_pred_indices so the order will be 0,1,...,n

    y_pred_indices = [index for sublist in [fold['teIndices'] for fold in xp_indices] for index in sublist]

    y_pred = list()

    for fold in xrange(len(xp_indices)) :
        with open(y_pred_fname%fold, 'r') as f_pred:
            content = f_pred.read().split()
            y_pred.extend(float(y) for y in content)
     
    y_pred_sorted = [y_pred[i] for i in y_pred_indices]

    # compute rmse using metrics :
    rmse = sklearn.metrics.mean_squared_error(y_true, y_pred_sorted)

    # output :
    with open(output_fname, 'a') as f_out:
        f_out.write("%d\n" %rmse)


def compute_ppv_sensitivity(causal_fname, selected_list, num_features):
    """ Compute PPV (Positive Predicted Values) = Accuracy = Precision
    and sensitivity (true positive rate) for all tasks.

    Arguments
    ---------
    causal_fname: filename
        File containing causal features (one line per task, space-separated).
    selected_list: list of lists
        List of lists of selected features (one list per task).
    num_features : int
        Total number of features

    Returns
    -------
    ppv_list: list
        List of Positive Predicted Values (PPV), task per task.
    tpr_list: list
        List of sensitivities (TPR), task per task.
    """

    ppv_list = []
    tpr_list = []

    # For each task, at the beginning, we consider that the features are neither causal nor predicted as such.
    # Then we change the state/status of the causal ones (these are y_true True),
    # and of those that have been predicted as such (these are y_pred True).
    # and we compute ppv and tpr based on these 2 sets.
    with open(causal_fname, 'r') as f:
        for line_idx, line in enumerate(f):

            y_true = [False]*num_features
            y_pred = [False]*num_features

            print 'line idx = ', line_idx
            y_true_indx_list = map(int, line.split())
            y_pred_indx_list = selected_list[line_idx]

            for y_true_indx in y_true_indx_list :
                y_true[y_true_indx] = True
            for y_pred_indx in y_pred_indx_list :
                y_pred[y_pred_indx] = True

            ppv_list.append( sklearn.metrics.accuracy_score(y_true, y_pred) )

            count_tpr = 0
            for i, j in zip(y_pred, y_true):
                if (i == j):
                    count_tpr += 1
            tpr_list.append( count_tpr / num_features)
    
    return ppv_list, tpr_list

    
class Framework(object):
    """ Setting up evaluation framework.

    Attributes
    ----------
    self.num_samples: int
        Number of samples.
    self.num_folds: int
        Number of cross-validation folds
    self.num_subsamples: int
        Number of subsamples (to evaluate stability)
    self.xp_indices: list of dictionaries
        fold_idx
        {
            'trIndices': list of train indices,
            'teIndices': list of test indices,
            'ssIndices': list of list of subsample indices
        }

    """
    def __init__(self, num_samples, num_folds, num_subsamples):
        """
        Parameters
        ----------
        num_samples: int
            Number of samples.
        num_folds: int
            Number of cross-validation folds
        num_subsamples: int
            Number of subsamples (to evaluate stability)
        """
        self.num_samples = num_samples
        self.num_folds = num_folds
        self.num_subsamples = num_subsamples
        self.xp_indices = [{'trIndices': list(), 'teIndices':list(), 'ssIndices':list()} for fold in xrange(num_folds)]
        
    def compute_indices(self, seed=None):
        """ Compute the cross-validation folds and subsample indices.

        Parameters
        ----------
        seed: {int, None}, optional
            random seed.
            Will always return the same with the same random seed.

        Modified attributes
        -------------------
        xp_indices: list of dictionaries
            fold_idx
            {
                'trIndices': list of train indices,
                'teIndices': list of test indices,
                'ssIndices': list of list of subsample indices
            }
        """

        # use sklearn.cross_validation
        
        # Generate cross-validation indices
        kf = cv.KFold(self.num_samples, n_folds=self.num_folds)
        for i, (train_index, test_index) in enumerate(kf):
            self.xp_indices[i]['trIndices'] = train_index.tolist()
            self.xp_indices[i]['teIndices'] = test_index.tolist()

            # For each train set, generate self.num_subsamples subsample sets of indices
            ss = cv.KFold(n=self.num_samples, n_folds=self.num_folds, shuffle=True, random_state=seed)
            for train_index, test_index in ss:
                self.xp_indices[i]['ssIndices'].append(train_index.tolist())
        
    def save_indices(self, data_dir, simu_id):
        """ Save the cross-validation folds and subsample indices to files.

        Parameters
        ----------

        Generated files
        ---------------
        For each fold_idx:
            <data_dir>/<simu_id>.<fold_idx>.trIndices:
                Space-separated list of training indices.
            <data_dir>/<simu_id>.<fold_idx>.teIndices:
                Space-separated list of test indices.
            For each subsample_idx:
                <data_dir>/<simu_id>.<fold_idx>.<ss_idx>.ssIndices
                    Space-separated lists of subsample indices,
                    one line per list / subsample.
        """
        # use np.savetxt 