"""evaluation_framework.py -- All that is needed to evaluate feature selection algorithms."""

import numpy as np
import sklearn
import sklearn.linear_model as lm
import sklearn.cross_validation as cv
import tables as tb
import subprocess
import shlex
import math


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
    cidx = -1.
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
    cidx = 0.
    for k1, sel1 in enumerate(sel_list[:-1]):
        # sel_list[:-1] to not take into account the last list.
        # avoid a problem with sel_list[k1+1:] when k1 is the last element,
        # that give an empty list overwise
        # the work is done at the second to last element anyway
        for sel2 in sel_list[k1+1:]:
            cidx += consistency_index(set(sel1), set(sel2), num_features)
    cidx = 2. / (len(sel_list) * (len(sel_list) - 1)) * cidx
    return cidx


def consistency_index_task(selection_fname, num_folds, num_tasks, num_features):
    """ Compute consistency indices between the features selected for each fold at each task

    Arguments
    ---------
    selection_fname : filename
        Template of path where were write list of selected features
    num_folds: int
        Total number of fold
    num_tasks: int
        Total number of task
    num_features: int
        Total number of features

    Return
    -------
    ci_list
        List of consistency indices between the features selected for each fold, task per task.
    """
    # In the curret repeat repo, there is a file of selected feature for each fold
    # in which each line is a task 
    # on each line, there is the space separated list of selected feature.
    # we want consistency indices between the features selected for each fold, task per task
    # so for each line in these files, we compute the ci between the features selected for each fold


    # As there are ~10 folds, there are 10 files.
    # they don't take a lot of memory
    # so it is ok to open them all :
    fold_f_list = []
    for fold_idx in xrange (num_folds) :
        f_sel = open(selection_fname %fold_idx, 'r')
        fold_f_list.append(f_sel)

    
    ci_list = []
    # For each task : 
    for task_idx in xrange (num_tasks):
        sel_list = []
        for f_sel in fold_f_list :
            # increment aline in each file
            content = f_sel.readline().split()
            # append lines content in sel_list
            sel_list.append(content)
        # compute the ci between the features selected for each fold at this current task
        ci =  consistency_index_k(sel_list, num_features)
        ci_list.append(ci)

    for f in fold_f_list :
        f.close()

    return ci_list


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
    #                    lbd, eta, 0, covariance_fname)
    # tt = sfan_solver.create_dimacs()
    # sfan_solver.run_maxflow()

    # But because cython output to screen is NOT caught by sys.stdout, 
    # we need to run this externally
    argum = ['/usr/bin/time', '--format=%M', 
             'python', 'multitask_sfan.py',
             '--num_tasks', str(num_tasks),
             '--networks', network_fname,
             '--node_weights']
    argum.extend(weights_fnames)
    argum.extend(params.split())
    argum.extend(['-m', '0'])
    print '+++'
    print argum
    p = subprocess.Popen(argum, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stdout=subprocess.PIPE -> something should read the output while the process is still running
    # stderr=subprocess.STDOUT : To also capture standard error in the result

    p_com = p.communicate()
    p_out = p_com[0].split("\n")
    p_err = p_com[1].split("\n")
    print p_com
    # Process the output to get lists of selected features
    sel_list = [[(int(x)-1) for x in line.split()] for line in p_out[2:2+num_tasks]]

    if not sel_list :
        #TODO : fix no sel_list issue#1
        #import pdb ; pdb.set_trace()
        print "WARNING : returned sel_list empty !! algo = st ; param = ", params
        sel_list = [[] for i in xrange(num_tasks)]

    # Process the standart output to get timing info 
    timing = '\n'.join(p_out[2+num_tasks:])

    # Process the standart error to get maxRSS info : 
    maxRSS = p_err[-2]

    return sel_list, timing, maxRSS
                 

def run_msfan_nocorr(num_tasks, network_fname, weights_fnames, params):
    """ Run multitask sfan (no precision/covariance matrix).

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
    argum = ['/usr/bin/time', '-f', '%M',
             'python', 'multitask_sfan.py',
             '--num_tasks', str(num_tasks),
             '--networks', network_fname,
             '--node_weights']
    argum.extend(weights_fnames)
    argum.extend(params.split())
    print '+++'
    print argum
    p = subprocess.Popen(argum, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    p_com = p.communicate()
    p_out = p_com[0].split("\n")
    p_err = p_com[1].split("\n")
    print p_com
    # Process the output to get lists of selected features
    
    sel_list = [[(int(x)-1) for x in line.split()] for line in p_out[3:3+num_tasks]]

    if not sel_list :
        #TODO : fix no sel_list issue#1
        #import pdb ; pdb.set_trace()
        print "WARNING : returned sel_list empty !! algo = np ; param = ", params
        sel_list = [[] for i in xrange(num_tasks)]

    # Process the output to get timing info 
    timing = '\n'.join(p_out[3+num_tasks:])

    # Process the outut to get maxRSS info : 
    maxRSS = p_err[-2]

    return sel_list, timing, maxRSS
                 

def run_msfan(num_tasks, network_fname, weights_fnames, covariance_fname, params):
    """ Run multitask sfan.

    Arguments
    ---------
    num_tasks: int
        Number of tasks. 
    network_fname: filename
        Path to the network file.
    weights_fnames: list of filenames
        List of paths to the network nodes files (one per task).
    covariance_fname: filename
        Path to the matrix of covariance (similarity) of tasks.
    params: string
        Hyperparameters, in the '-l <lambda> -e <eta> -m <mu>' format.

    Returns
    -------
    sel_list: list of lists
        For each task, a list of selected features, as indices,
        STARTING AT 0.
    """
    argum = ['/usr/bin/time', '-f', '%M',
             'python', 'multitask_sfan.py',
             '--num_tasks', str(num_tasks),
             '--networks', network_fname,
             '--node_weights']
    argum.extend(weights_fnames)
    argum.extend(['--covariance_matrix', covariance_fname])
    argum.extend(params.split())
    print '+++'
    print argum
    p = subprocess.Popen(argum, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    p_com = p.communicate()
    p_out = p_com[0].split("\n")
    p_err = p_com[1].split("\n")
    print p_com
    # Process the output to get lists of selected features
    sel_list = [[(int(x)-1) for x in line.split()] for line in p_out[3:3+num_tasks]]

    if not sel_list : 
        #TODO : fix no sel_list issue#1
        #import pdb ; pdb.set_trace() 
        print "WARNING : returned sel_list empty !! algo = msfan ; param = ", params
        sel_list = [[] for i in xrange(num_tasks)]

    # Process the output to get timing info 
    timing = '\n'.join(p_out[3+num_tasks:])

    # Process the outut to get maxRSS info : 
    maxRSS = p_err[-2]

    return sel_list, timing, maxRSS
                 

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
        Optimal parameters, leading to highest consistency index mean
        of features selected for each subsample for each task
        => params leading to the best ci mean.
    """
    opt_params = ''
    opt_ci_mean = -1 # set to -1 because it is the worst case ci value 
    for (params, selected_dict_p) in selected_dict.iteritems():
        ci_list = [] #list of ci, one per task, computed with current params
        for (task_idx, sel_list) in selected_dict_p.iteritems():
            ci_of_current_task = consistency_index_k(sel_list, num_features)
            ci_list.append(ci_of_current_task) 
        ci_mean = np.mean(ci_list)
        if ci_mean >= opt_ci_mean:
            opt_ci_mean = ci_mean
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

    #----------------------------------------
    # Read data : 

    if not selected_features :
        # Safeguard for when SFAN returns empty list
        # Avoid not allowed empty selections 
        #TODO : fix no sel_list issue#1
        #import pdb; pdb.set_trace() 
        print ('WARNING : no features was selected on this fold -> give NA predictions') 
        preds = np.array([np.nan ] * len(te_indices) )
    else :
        # read genotypes : 
        with tb.open_file(genotype_fname, 'r') as h5f:
            table = h5f.root.Xtr
            # table.shape : 
            # For each feature (in line), 
            # there is the genotype of each sample (in column)
            X = table[selected_features, :]

        Xtr = [X[:,tr] for tr in tr_indices]
        Xte = [X[:,te] for te in te_indices]

        # read phenotypes : 
        with open(phenotype_fname, 'r') as f:
            # continuous phenotype for each sample (in line)
            y = f.read().split()
            y = [float(item) for item in y]

        ytr = [ y[tr] for tr in tr_indices]
        #----------------------------------------

        # Instantiate a ridge regression
        model = lm.RidgeCV()

        # Train the ridge regression on the training set
        model.fit(Xtr, ytr)

        #----------------------------------------

        # Make predictions on the test set
        preds = model.predict(Xte)

    # Save predictions
    np.savetxt(output_fname, preds, fmt='%.3e')




def compute_ridge_selected_RMSE(phenotype_fnames, y_pred_template, xp_indices):
    """ Compute RMSE (Root Mean Squared Error)

    Arguments
    ---------
    phenotype_fnames: aray of filename
        Path to phenotype datas.
        (one path per task)
    y_pred_template: string
        Template of path where were write list of predictions on the test set
    xp_indices: list of dictionaries
        fold_idx
        {
            'trIndices': list of train indices,
            'teIndices': list of test indices,
            'ssIndices': list of list of subsample indices
        }
    

    Return
    -------
    rmse_list:
        List of rmse task per task.
    """
    rmse_list = []
    for task_idx, phenotype_fname in enumerate( phenotype_fnames ) :
        #print "\n\n\n\n==== tache num %d" %task_idx
        # For n inds :
        # RMSE = sqrt { (1/n)  [sum from m=1 to n : (ypred_m - ytrue_m)^2 ]  }

        # read all_y_true :
        #print '\ni read phenotype_fnames[task_idx = %d] = %s' %(task_idx, phenotype_fnames[task_idx]) 
        with open(phenotype_fname, 'r') as f_true:
            all_y_true = [float(y) for y in f_true.read().split()]
        #print "\nall_y_true = "
        #print all_y_true
        # read all_y_pred :
        # predictions were made one by one, in order : [fold['teIndices'] for fold in xp_indices]
        # we open each file (one per fold) and append predicted phenotypes
        # then when sort them using all_y_pred_indices so the order will be 0,1,...,n

        all_y_pred_indices = [index for sublist in [fold['teIndices'] for fold in xp_indices] for index in sublist]
        all_y_pred = list()

        num_folds = len(xp_indices)
        for fold_idx in xrange(num_folds) :
            with open(y_pred_template%(fold_idx, task_idx), 'r') as f_pred:
                content = f_pred.read().split()
                all_y_pred.extend(float(y) for y in content)
         
        all_y_pred_sorted = [all_y_pred[i] for i in all_y_pred_indices]
        #print "\n all_y_pred_sorted = "
        #print all_y_pred_sorted
        # compute rmse using metrics : 
        # wanted to use : rmse = sklearn.metrics.mean_squared_error(all_y_true, all_y_pred_sorted)
        # be if all_y_pred_sorted have NaN, there is a problem
        # -> compute RMSE without NaN ind

        #TODO : is it possible to have NA only for some ind ? 
        #if not : use if np.isnan(all_y_pred_sorted).any() 
        #and don't compute RMSE for this task without using not_NaN_idx
        not_NaN_idx = np.where(~ np.isnan(all_y_pred_sorted) )[0]
        if not not_NaN_idx.size : # if not_NaN_idx empty -> if there is only NaNs in all_y_preds
            rmse = np.NaN
        else : 
            not_NaN_y_true =  [all_y_true[i] for i in not_NaN_idx]
            not_NaN_y_pred_sorted = [all_y_pred_sorted[i] for i in not_NaN_idx]
            rmse = math.sqrt (sklearn.metrics.mean_squared_error(not_NaN_y_true, not_NaN_y_pred_sorted) )

        #print "rmse = %f" % rmse
        rmse_list.append(rmse)

    # return :
    return rmse_list


def evaluate_classification(causal_features, selected_features, num_features):
    """ Compute metrics scoring classification, for all tasks.

    Arguments
    ---------
    causal_features:  list of lists
        List of lists of real causal features (one list per task).
    selected_features: list of lists
        List of lists of selected features (one list per task).
    num_features : int
        Total number of features

    Returns
    -------
    acc_list: list
        List of Accuracy, task per task.
        fraction of correct predictions = predictions matching observations.
    mcc_list: list
        List of Matthews correlation coefficient task per task.
        Better than Accuracy because classes are of very different sizes.
    pre_list: list
        List of Positive Precision = Predicted Values (PPV), task per task.
        = TP / (TP + FP)
    spe_list: list
        List of Recall = Sensitivity =  True Positive Rate (TPR), task per task.
        = TP / (TP + FN)
    """
    acc_list = []
    mcc_list = []
    pre_list = []
    spe_list = []
    
    
    # For each task, 
    for task_idx in xrange(len(causal_features)):

        # at the beginning, we consider that the features are 
        # neither causal...
        y_true = [False]*num_features
        # ... nor predicted as such.
        y_pred = [False]*num_features

        # Then we change the status of the causal ones 
        # (these are y_true True),
        for y_true_idx in causal_features[task_idx] :
            y_true[y_true_idx] = True
        # and of those that have been predicted as such 
        # (these are y_pred True). 
        for y_pred_idx in selected_features[task_idx] :
            y_pred[y_pred_idx] = True

        # and we compute 
        #   - accuracy_score
        #   - matthews_corrcoef
        #   - precision_score
        #   - recall_score
        # based on these 2 sets : 
        acc_list.append( sklearn.metrics.accuracy_score    (y_true, y_pred) )
        mcc_list.append( sklearn.metrics.matthews_corrcoef (y_true, y_pred) )
        pre_list.append( sklearn.metrics.precision_score   (y_true, y_pred) )
        spe_list.append( sklearn.metrics.recall_score      (y_true, y_pred) )

    return acc_list, mcc_list, pre_list, spe_list


    
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
        kf = cv.KFold(self.num_samples, n_folds=self.num_folds, shuffle = True, random_state=seed)
        for fold_idx, (train_indices_f, test_indices_f) in enumerate(kf):
            #print fold_idx, train_indices_f, test_indices_f
            # Generate cross-validation indices
            self.xp_indices[fold_idx]['trIndices'] = train_indices_f.tolist()
            self.xp_indices[fold_idx]['teIndices'] = test_indices_f.tolist()
            # For each train set, generate self.num_subsamples subsample sets of indices (90% of the train_set_f)
            for i_ss in xrange(self.num_subsamples) : 
                train_indices_ss, test_indices_ss = cv.train_test_split(train_indices_f, train_size=0.9)
                self.xp_indices[fold_idx]['ssIndices'].append( train_indices_ss.tolist() ) 

        
    def save_indices(self, out_dir, simu_id):
        """ Save the cross-validation folds and subsample indices to files.

        Parameters
        ----------
        out_dir : dir path
            fold where indices have to be saved
        simu_id :  string
            Name of the simulation, to be used to name files.
        
        Generated files
        ---------------
        For each fold_idx:
            <out_dir>/<simu_id>.fold<fold_idx>.trIndices:
                Space-separated list of training indices.
            <out_dir>/<simu_id>.fold<fold_idx>.teIndices:
                Space-separated list of test indices.
            For each subsample_idx:
                <out_dir>/<simu_id>.fold<fold_idx>.ss<ss_idx>.ssIndices
                    Space-separated lists of subsample indices,
                    one line per list / subsample.
        """
        # use np.savetxt ??? why ? 
        trIndices_fname = out_dir+'/'+simu_id+'.fold%d.trIndices'
        teIndices_fname = out_dir+'/'+simu_id+'.fold%d.teIndices'
        ssIndices_fname = out_dir+'/'+simu_id+'.fold%d.ss%d.ssIndices'
        for fold_idx in xrange(self.num_folds) : 
            with open(trIndices_fname %(fold_idx), 'w') as trIndices_f : 
                trIndices_f.write(  " ".join(str(i) for i in self.xp_indices[fold_idx]["trIndices"] ) )
            #np.savetxt(trIndices_fname %(fold_idx), self.xp_indices[fold_idx]["trIndices"], delimiter=' ', fmt='%d')
            with open(teIndices_fname %(fold_idx),'w') as teIndices_f : 
                teIndices_f.write(  " ".join(str(i) for i in self.xp_indices[fold_idx]["teIndices"] ) )
            #np.savetxt(teIndices_fname %(fold_idx), self.xp_indices[fold_idx]["teIndices"], delimiter=' ', fmt='%d')
            for ss_idx in xrange(self.num_subsamples) :
                with open(ssIndices_fname  %(fold_idx,ss_idx), 'w') as ssIndices_f:
                    ssIndices_f.write(  " ".join(str(i) for i in self.xp_indices[fold_idx]["ssIndices"][ss_idx] ) )
                #np.savetxt(ssIndices_fname  %(fold_idx,ss_idx), self.xp_indices[fold_idx]["ssIndices"][ss_idx], delimiter=' ', fmt='%d')
