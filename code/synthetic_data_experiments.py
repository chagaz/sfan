# -*- coding: utf-8 -*-
"""synthetic_data_experiments.py -- Run validation experiments on synthetic data

In this version all experiments are run sequentially.
"""

DEBUG_MODE = True
DATA_GEN = False # have to gene dat or not ?
SEQ_MODE = False

# Importing local libraries first,
# because otherwise Error in `python': free(): invalid pointer
import multitask_sfan
import evaluation_framework as ef
import generate_data
import plot

import argparse
import logging
import os
import numpy as np
import scipy.stats as st
import subprocess
import sys
import tables as tb
import tempfile
import shutil
import shlex
import glob


def get_arguments_values(): 
    """ Use argparse module to get arguments values.

    Returns
    -------
    args : Namespace object
        Namespace object populated with arguments converted from strings to objects 
        and assigned as attributes of the namespace object.
        They store arguments values (str or int, according to the specifications in
        the code)
    """
    help_str = "Validation experiments on synthetic data"
    parser = argparse.ArgumentParser(description=help_str,add_help=True)
    parser.add_argument("-k", "--num_tasks", help="Number of tasks", type=int)
    parser.add_argument("-m", "--num_features", help="Number of features",
                        type=int)
    parser.add_argument("-n", "--num_samples", help="Number of samples",
                        type=int)
    parser.add_argument("-r", "--num_repeats", help="Number of repeats",
                        type=int)
    parser.add_argument("-f", "--num_folds", help="Number of CV folds",
                        type=int)
    parser.add_argument("-s", "--num_subsamples", help="Number of subsamples",
                        type=int)
    parser.add_argument("data_dir", help="Simulated data directory")
    parser.add_argument("resu_dir", help="Results directory")
    parser.add_argument("simu_id", help="Simulation name")
    parser.add_argument("-v", "--verbose", help="Turn on detailed info log",
                        action='store_true')
    return parser.parse_args()

def check_arguments_integrity(args): 
    """ Check integrity of arguments pass through the command line. 

    Parameters
    ----------
    args : namespace object
        Its attributes are arguments names 
        and contain arguments values. 

    Side effects
    ------------
    Exit the script printing why if one arg is not integrous. 
    """
    try:
        assert(args.num_tasks >= 1)
    except AssertionError:
        logging.error("There must be at least one task specified.\n")
        logging.error("Use --help for help.\n")
        sys.exit(-1)
    
    try:
        assert(args.num_features >= generate_data.NUM_CAUSAL_TOTAL)
    except AssertionError:
        logging.error("The number of features must be larger than" + \
                      " NUM_CAUSAL_TOTAL (%d).\n" % \
                      generate_data.NUM_CAUSAL_TOTAL)
        logging.error("Use --help for help.\n")
        sys.exit(-1)

    try:
        assert(args.num_samples > 0)
    except AssertionError:
        logging.error("The number of samples must be strictly positive\n")
        logging.error("Use --help for help.\n")
        sys.exit(-1)

    try:
        assert(args.num_repeats > 0)
    except AssertionError:
        logging.error("The number of repeats must be strictly positive\n")
        logging.error("Use --help for help.\n")
        sys.exit(-1)

    try:
        assert(args.num_folds > 0)
    except AssertionError:
        logging.error("The number of cross-validation folds must be strictly positive\n")
        logging.error("Use --help for help.\n")
        sys.exit(-1)

    try:
        assert(args.num_subsamples > 0)
    except AssertionError:
        logging.error("The number of subsamples must be strictly positive\n")
        logging.error("Use --help for help.\n")
        sys.exit(-1)

    # Verbose
    if args.verbose:
        logging.basicConfig(format="[%(levelname)s] %(message)s",
                            level=logging.DEBUG)
        logging.info("Verbose output.")
    else:
        logging.basicConfig(format="%[(levelname)s] %(message)s")

def get_integrous_arguments_values(): 
    """ Get arguments passed through command line and check their integrity.
    Returns
    -------
    args : namespace object
        Its attributes are arguments names 
        and contain arguments values (str or int according to code specifications). 
    """
    args = get_arguments_values()
    check_arguments_integrity(args)
    return args

def create_dir_if_not_exists(dir_name): 
    """ Create a dir named <dir_name> if it does not already exists. 

    Parameters
    ----------
    dir_name : str
        Path of the wanted new dir. 

    Side effects
    ------------
    Create a dir named <dir_name> if it does not already exists.
    """
    if not os.path.isdir(dir_name):
        logging.info("Creating %s\n" % dir_name)
        try: 
            os.makedirs(dir_name)
        except OSError:
            if not os.path.isdir(dir_name):
                raise

def get_analysis_files_names(resu_dir, simu_id): 
    """ Give analysis files names. 

    Parameters
    ----------
    resu_dir: filename
        Path of the directory in which to save the results.
    simu_id: string
        Name of the simulation, to be used to name files.

    Returns
    -------
    analysis_files : dictionary
        key : <measure>_<algo> 
        value : filename
    """
    # - to hold PPV values
    ppv_st_fname = '%s/%s.sfan.ppv' % (resu_dir, simu_id)
    ppv_np_fname = '%s/%s.msfan_np.ppv' % (resu_dir, simu_id)
    ppv_fname = '%s/%s.msfan.ppv' % (resu_dir, simu_id)
    # - to hold sensitivity = TPR values
    tpr_st_fname = '%s/%s.sfan.sensitivity' % (resu_dir, simu_id)    
    tpr_np_fname = '%s/%s.msfan_np.sensitivity' % (resu_dir, simu_id)
    tpr_fname = '%s/%s.msfan.sensitivity' % (resu_dir, simu_id)
    # - to hold consistency values
    ci_st_fname = '%s/%s.sfan.consistency' % (resu_dir, simu_id)
    ci_np_fname = '%s/%s.msfan_np.consistency' % (resu_dir, simu_id)
    ci_fname = '%s/%s.msfan.consistency' % (resu_dir, simu_id)
    # - to hold RMSE values
    rmse_st_fname = '%s/%s.sfan.rmse' % (resu_dir, simu_id)
    rmse_np_fname = '%s/%s.msfan_np.rmse' % (resu_dir, simu_id)
    rmse_fname = '%s/%s.msfan.rmse' % (resu_dir, simu_id)
    # - to hold timing values
    timing_st_fname = '%s/%s.sfan.timing' % (resu_dir, simu_id)
    timing_np_fname = '%s/%s.msfan_np.timing' % (resu_dir, simu_id)
    timing_fname = '%s/%s.msfan.timing' % (resu_dir, simu_id)
    # - to hold maxRSS values 
    maxRSS_st_fname = '%s/%s.sfan.maxRSS' % (resu_dir, simu_id)
    maxRSS_np_fname = '%s/%s.msfan_np.maxRSS' % (resu_dir, simu_id)
    maxRSS_fname = '%s/%s.msfan.maxRSS' % (resu_dir, simu_id)
    
    analysis_files = {
        'ppv_st':ppv_st_fname,
        'ppv_msfan_np':ppv_np_fname,
        'ppv_msfan':ppv_fname,
        'tpr_st':tpr_st_fname,
        'tpr_msfan_np':tpr_np_fname,
        'tpr_msfan':tpr_fname,
        'ci_st':ci_st_fname ,
        'ci_msfan_np':ci_np_fname ,
        'ci_msfan':ci_fname, 
        'rmse_st':rmse_st_fname ,
        'rmse_msfan_np':rmse_np_fname ,
        'rmse_msfan':rmse_fname,
        'timing_st':timing_st_fname,
        'timing_msfan_np':timing_np_fname,
        'timing_msfan':timing_fname,
        'maxRSS_st':maxRSS_st_fname,
        'maxRSS_msfan_np':maxRSS_np_fname,
        'maxRSS_msfan':maxRSS_fname
    }

    return analysis_files

def determine_hyperparamaters(genotype_fname, phenotype_fnames, network_fname, precision_fname, args):
    """ Determine hyperparameters. 

    Parameters
    ----------

    genotype_fname: filename
        Path to genotype data.
    phenotype_fname: filename
        Path to phenotype data.
    network_fname: filename
        Path to the network file.
    precision_fname: filename
        Path to the precision / correlation matrix file. 
    args : Namespace object
        Its attributes are arguments names 
        and contain arguments values (str or int according to code specifications).
    Returns
    -------
    lbd_eta_values: list of strings
        Values of hyperparameters, in the format:
        "-l <lambda -e <eta>". -> for Sfan
    lbd_eta_mu_values_np: list of strings
        Values of hyperparameters, in the format:
        "-l <lambda -e <eta> -m <mu>". -> for MSfan_np 
    lbd_eta_mu_values: list of strings
        Values of hyperparameters, in the format:
        "-l <lambda -e <eta> -m <mu>". -> for MSfan

    """
    # Define the grid of hyperparameters
    # see paper/tech_note
    # Randomly sample 50% of the data
    tmp_scores_f_list = []
    with tb.open_file(genotype_fname, 'r') as h5f:
        Xtr = h5f.root.Xtr[:, :]

        # Define subsample of 50% of the data
        sample_indices = range(args.num_samples)
        np.random.shuffle(sample_indices)
        sample_indices = sample_indices[:(args.num_samples/2)]
        Xtr = Xtr[:, sample_indices]

        # Compute scores for the subsample
        for task_idx in range(args.num_tasks):
            # Read phenotype
            y = np.loadtxt(phenotype_fnames[task_idx])[sample_indices]

            # Compute feature-phenotype correlations
            r2 = [st.pearsonr(Xtr[feat_idx, :].transpose(), y)[0]**2 \
                  for feat_idx in range(args.num_features)]

            # Create temporary file of name tmp_fname 
            fd, tmp_fname = tempfile.mkstemp()

            # Save to temporary file
            np.savetxt(tmp_fname, r2, fmt='%.3e')

            # Append temporary file to list
            tmp_scores_f_list.append(tmp_fname)

    # Compute grid (WARNING: STILL NOT WORKING WELL)
    sfan_ = multitask_sfan.Sfan(args.num_tasks, [network_fname],
                                tmp_scores_f_list, 1, 1, 1,
                                precision_matrix_f=precision_fname)
    # not 0, 0, 0 but 1, 1, 1 (or anything else > 0)
    # because if mu = 0, the matrix is not use 
    # -> no matrix, no phi, etc. 

    # for Sfan & MSfan np : 
    lbd_eta_mu_values_np = sfan_.compute_hyperparameters_range(num_values=5)
    lbd_eta_values = [" ".join(plist.split()[:-2]) \
                      for plist in lbd_eta_mu_values_np]
    # for MSfan : 
    lbd_eta_mu_values = sfan_.compute_hyperparameters_range_multiscones(num_values=5)

    # Delete temporary files from tmp_scores_f_list
    for fname in tmp_scores_f_list:
        os.remove(fname)

    return lbd_eta_values, lbd_eta_mu_values_np, lbd_eta_mu_values


def run_fold(fold_idx, args, lbd_eta_values, lbd_eta_mu_values, indices, genotype_fname, network_fname , precision_fname , causal_fname, phenotype_fnames, scores_fnames, resu_dir):
    """TODO
    """

    if not DEBUG_MODE : 
        logging.info ("======== Feature selection :")
        #XXX DEBUG ???

        # Inititalize dictionary to store selected features
        # sf_dict is a nested dictionary, indexed by
        #   - value of the parameters
        #   - value of task_idx
        #   - subsample idx
        # sf is a list of lists of selected features
        # (one per subsample iteration)
        # you get a specific sf from sf_dict by querying
        # sf_dict[params][task_idx]
        sf_st_dict = {}       # single task
        sf_np_dict = {}       # not using precision matrix
        sf_dict = {}          # using precision matrix
        for params in lbd_eta_values:
            sf_st_dict[params] = {}
            for task_idx in range(args.num_tasks):
                sf_st_dict[params][task_idx] = []
        for params in lbd_eta_mu_values_np:
            sf_np_dict[params] = {}
            for task_idx in range(args.num_tasks):
                sf_np_dict[params][task_idx] = []
        for params in lbd_eta_mu_values:
            sf_dict[params] = {}
            for task_idx in range(args.num_tasks):
                sf_dict[params][task_idx] = []

        #process_time files template : 
        process_time_file_template = resu_dir+'/'+args.simu_id+'.%s.fold_'+str(fold_idx)+'ss.process_time'
        #max RSS files template : 
        max_RSS_file_template = resu_dir+'/'+args.simu_id+'.%s.fold_'+str(fold_idx)+'ss.max_RSS'

        for ss_idx in range(args.num_subsamples):
            logging.info ("========                        SS : %d" % ss_idx)
            # Get samples
            sample_indices = evalf.xp_indices[fold_idx]['ssIndices'][ss_idx]
        
            # Generate sample-specific network scores from phenotypes and genotypes
            tmp_weights_f_list = [] # to hold temp files storing these scores
            print ('DEBUG je vais lire un pytables file')
            with tb.open_file(genotype_fname, 'r') as h5f:
                print ("DEBUG j'ai ouvert un pytables file")
                Xtr = h5f.root.Xtr[:, sample_indices]
                for task_idx in range(args.num_tasks):
                    # Read phenotype
                    y = np.loadtxt(phenotype_fnames[task_idx])[sample_indices]
        
                    # Compute feature-phenotype correlations
                    r2 = [st.pearsonr(Xtr[feat_idx, :].transpose(), y)[0]**2 \
                          for feat_idx in range(args.num_features)]
        
                    # Save to temporary file tmp_weights_f_list[task_idx]
                    # Create temporary file of name tmp_fname (use tempfile)
                    fd, tmp_fname = tempfile.mkstemp()
                    # Save to temporary file
                    np.savetxt(tmp_fname, r2, fmt='%.3e')
                    # Append temporary file to list
                    tmp_weights_f_list.append(tmp_fname)
        
            for params in lbd_eta_values:
                logging.info("========                        lbd_eta_values : "+ `params`)
                # Select features with single-task sfan
                logging.info("                                   run_sfan")
                sel_, timing, max_RSS = ef.run_sfan(args.num_tasks, network_fname,
                                   tmp_weights_f_list, params)
                if not sel_ : import pdb; pdb.set_trace() #DEBUG
                # Store selected features in the dictionary
                for task_idx, sel_list in enumerate(sel_):
                    sf_st_dict[params][task_idx].append(sel_list)
                #Store process time
                process_time = timing.split()[-1]
                fname= process_time_file_template % 'sfan'
                with open(fname, 'a') as f:
                    f.write("%s\n" % process_time)
                #Store max RSS : 
                fname= max_RSS_file_template % 'sfan'
                with open(fname, 'a') as f:
                    f.write("%s\n" % max_RSS)
            
            for params in lbd_eta_mu_values_np:
                logging.info("========                        lbd_eta_mu_values_np"+ `params`)
                # Select features with multi-task (no precision) sfan
                logging.info("                                   run_msfan_nocorr")
                sel_ , timing, max_RSS = ef.run_msfan_nocorr(args.num_tasks, network_fname,
                                           tmp_weights_f_list, params)
                if not sel_ : import pdb; pdb.set_trace()#DEBUG
                # Store selected features in the dictionary
                for task_idx, sel_list in enumerate(sel_):
                    sf_np_dict[params][task_idx].append(sel_list)
                #Store process time
                process_time = timing.split()[-1]
                fname= process_time_file_template % 'msfan_np'
                with open(fname, 'a') as f:
                    f.write("%s\n" % process_time)
                #Store max RSS : 
                fname= max_RSS_file_template %  'msfan_np'
                with open(fname, 'a') as f:
                    f.write("%s\n" % max_RSS)

            for params in lbd_eta_mu_values:
                logging.info("========                        lbd_eta_mu_values"+ `params`)
                # Select features with multi-task sfan
                logging.info("                                   run_msfan")
                sel_, timing, max_RSS = ef.run_msfan(args.num_tasks, network_fname,
                                    tmp_weights_f_list, precision_fname,
                                    params)
                if not sel_ : import pdb; pdb.set_trace() #DEBUG                                      
                # Store selected features in the dictionary
                for task_idx, sel_list in enumerate(sel_):
                    sf_dict[params][task_idx].append(sel_list)
                #Store process time
                process_time = timing.split()[-1]
                fname= process_time_file_template % 'msfan'
                with open(fname, 'a') as f:
                    f.write("%s\n" % process_time)
                #Store max RSS : 
                fname= max_RSS_file_template % 'msfan'
                with open(fname, 'a') as f:
                    f.write("%s\n" % max_RSS)


          
            # Delete the temporary files stored in tmp_weights_f_list
            for fname in tmp_weights_f_list:
                os.remove(fname)
        
        # END for ss_idx in range(args.num_subsamples)
        
        # XXX DEBUG : 
        sf_np_dict = sf_dict
        #-----------------------------------   
        # Get optimal parameter values for each algo.
        # ??? some lists are empty, is it normal ??? 
        logging.info( "======== Get opt params")
        opt_params_st = ef.get_optimal_parameters_from_dict(sf_st_dict, args.num_features)
        print 'opt param st ', opt_params_st
        opt_params_np = ef.get_optimal_parameters_from_dict(sf_np_dict, args.num_features)
        print 'opt param np ', opt_params_np
        opt_params = ef.get_optimal_parameters_from_dict(sf_dict, args.num_features)
        print 'opt params ', opt_params

    else : 
        # XXX DEBUG : ??? : 
        opt_params_st = "-l 5.97e-03 -e 1.60e-02"
        opt_params_np = "-l 2.94e-03 -e 7.65e-03 -m 1.42e-02"
        opt_params =    "-l 2.94e-03 -e 7.65e-03 -m 1.42e-02"

    # For each algorithm, save optimal parameters to file
    # Single task
    fname = '%s/%s.sfan.fold_%d.parameters' % (resu_dir, args.simu_id, fold_idx)
    with open(fname, 'w') as f:
        f.write(opt_params_st)
    # Multitask (no precision)
    fname = '%s/%s.msfan_np.fold_%d.parameters' % (resu_dir, args.simu_id, fold_idx)
    with open(fname, 'w') as f:
        f.write(opt_params_np)
    # Multitask (precision)
    fname = '%s/%s.msfan.fold_%d.parameters' % (resu_dir, args.simu_id, fold_idx)
    with open(fname, 'w') as f:
        f.write(opt_params)
    #------------------------------------------------------------------


    #------------------------------------------------------------------
    logging.info( "======== Features selection using all training set and opt param")
    #------
    # For each algorithm, run algorithms again to select features,
    # (got a list of list : list of selected features for each task)
    # using the whole training set (i.e. scores_fnames)
    # and optimal parameters.
    logging.info("          run st")
    selected_st, timing_st, maxRSS_st = ef.run_sfan(args.num_tasks, network_fname,
                               scores_fnames, opt_params_st)
    logging.info("          run np")
    selected_np, timing_np, maxRSS_np = ef.run_msfan_nocorr(args.num_tasks, network_fname,
                                       scores_fnames, opt_params_np)
    logging.info("          run msfan")
    selected, timing,maxRSS = ef.run_msfan(args.num_tasks, network_fname,
                                scores_fnames, precision_fname,
                                opt_params)

    #------
    # For each algorithm, save timing to file
    # Single task 
    with open(analysis_files['timing_st'], 'a') as f:
        f.write("%s\n" % timing_st)
    # Multitask (no precision)
    with open(analysis_files['timing_msfan_np'], 'a') as f:
        f.write("%s\n" % timing_np)
    # Multitask (precision)
    with open(analysis_files['timing_msfan'], 'a') as f:
        f.write("%s\n" % timing)


    #------
    # For each algorithm, save maxRSS to file
    # Single task 
    with open(analysis_files['maxRSS_st'], 'a') as f:
        f.write("%s\n" % maxRSS_st)
    # Multitask (no precision)
    with open(analysis_files['maxRSS_msfan_np'], 'a') as f:
        f.write("%s\n" % maxRSS_np)
    # Multitask (precision)
    with open(analysis_files['maxRSS_msfan'], 'a') as f:
        f.write("%s\n" % maxRSS)


    #------
    # For each algorithm, save selected features to file
    # Single task
    fname = '%s/%s.sfan.fold_%d.selected_features' % \
            (resu_dir, args.simu_id, fold_idx)
    
    with open(fname, 'w') as f:
        for selected_features_list in selected_st:
            f.write("%s\n" % ' '.join(str(x) for x in selected_features_list))
                                    # selected_features_list is a list of int 
                                    # that have to be cast as string so we can join them
    # Multitask (no precision)
    fname = '%s/%s.msfan_np.fold_%d.selected_features' % \
            (resu_dir, args.simu_id, fold_idx)
    
    with open(fname, 'w') as f:
        for selected_features_list in selected_np:
            f.write("%s\n" % ' '.join(str(x) for x in selected_features_list))
    # Multitask (precision)
    fname = '%s/%s.msfan.fold_%d.selected_features' % \
            (resu_dir, args.simu_id, fold_idx)
    
    with open(fname, 'w') as f:
        for selected_features_list in selected:
            f.write("%s\n" % ' '.join(str(x) for x in selected_features_list))
    
    #------
    # For each algorithm, and for each task, compute PPV
    # and sensitivity, and save to ppv_fname, tpr_fname
    #                                          v-- = fold_id
    # print in files : %s/repeat%d/%s.sfan.fold*.ppv" % (resu_dir, repeat_idx, simu_id)
    # with '%.2f' precision
    
    # Files structure : 
    # 1 line per repeat
    # on each line : valTask1, valTask2, ... valTaskn for each fold
    
    # For the current repeat and the current fold, 
    # ppv_list ant tpr_list and list of ppv ant tpr respectively
    # for each task



    ppv_template_f_name = str(resu_dir)+"/"+str(args.simu_id)+".%s.fold_"+str(fold_idx)+".ppv"
    tpr_template_f_name = str(resu_dir)+"/"+str(args.simu_id)+".%s.fold_"+str(fold_idx)+".tpr"

    # Single task
    ppv_list_st, tpr_list_st = ef.compute_ppv_sensitivity(causal_fname,
                                                    selected_st,
                                                    args.num_features)
    with open (ppv_template_f_name %"sfan", 'w') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in ppv_list_st]))
    with open (tpr_template_f_name %"sfan", 'w') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in ppv_list_st]))

    # Multitask (no precision)
    ppv_list_np, tpr_list_np = ef.compute_ppv_sensitivity(causal_fname,
                                                    selected_np,
                                                    args.num_features)
    with open (ppv_template_f_name %"msfan_np", 'w') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in ppv_list_np]))
    with open (tpr_template_f_name %"msfan_np", 'w') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in tpr_list_np]))

    # Multitask (precision)
    ppv_list_msfan, tpr_list_msfan = ef.compute_ppv_sensitivity(causal_fname,
                                                    selected,
                                                    args.num_features)
    with open (ppv_template_f_name %"msfan", 'w') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in ppv_list_msfan]))
    with open (tpr_template_f_name %"msfan", 'w') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in tpr_list_msfan]))

    #------------------------------------------------------------------


def run_repeat(repeat_idx, args, analysis_files):
    """ Run the repeat n° <repeat_idx>.

    Parameters
    ----------
    repeat_idx : int 
        Index of the current repeat. 
    args : Namespace object
        Its attributes are arguments names 
        and contain arguments values (str or int according to code specifications).
    analysis_files: dictionary
        key : <measure>_<algo> 
        value : filename

    """

    logging.info ("=============== REPETITION : %d" %repeat_idx)

    # Create <resu_dir>/repeat_<repeat_id> if it does not exist
    resu_dir = "%s/repeat_%d" % (args.resu_dir, repeat_idx)
    create_dir_if_not_exists(resu_dir)


    #-------------------------------------------------------------------------
    # Data generation : 
    logging.info( "======== Data generation")
    # Data generation : 

    # Instantiate data generator
    data_dir = '%s/repeat_%d' % (args.data_dir, repeat_idx)

    if DATA_GEN : 
        data_gen = generate_data.SyntheticDataGenerator(args.num_tasks,
                                                        args.num_features,
                                                        args.num_samples,
                                                        data_dir,
                                                        args.simu_id)
        # Generate modular data
        data_gen.generate_modular()

    # Name of data files
    # "Hard-coded here", but maybe edit generate_modular
    # to return the names of these files
    genotype_fname = '%s/%s.genotypes.txt' % (data_dir, args.simu_id)
    network_fname = '%s/%s.network.dimacs' % (data_dir, args.simu_id)
    precision_fname = '%s/%s.task_similarities.txt' % (data_dir,
                                                         args.simu_id)
    causal_fname = '%s/%s.causal_features.txt' % (data_dir, args.simu_id)
    phenotype_fnames = ['%s/%s.phenotype_%d.txt' % \
                        (data_dir, args.simu_id, task_idx) \
                        for task_idx in range(args.num_tasks)]
    scores_fnames = ['%s/%s.scores_%d.txt' % \
                     (data_dir, args.simu_id, task_idx) \
                     for task_idx in range(args.num_tasks)]

    
    #-------------------------------------------------------------------------


    #-------------------------------------------------------------------------
    logging.info("======== Ef setup")

    # Instantiate evaluation framework
    evalf = ef.Framework(args.num_samples, args.num_folds,
                         args.num_subsamples)

    # Compute cross-validation folds and subsample indices
    evalf.compute_indices()

    # Save cross-validation folds and subsample indices to file
    evalf.save_indices(data_dir, args.simu_id)


    #-------------------------------------------------------------------------
    # Looking for optimal parameters : 

    logging.info ("======== Defining grid of hyperparameters")
    lbd_eta_values, lbd_eta_mu_values_np, lbd_eta_mu_values  = determine_hyperparamaters( genotype_fname, 
                                                                    phenotype_fnames,
                                                                    network_fname,
                                                                    precision_fname,
                                                                    args)

    # and save them :
    hyperparam_fname = '%s/%s.hyperparameters.txt' % (data_dir, args.simu_id)
    with open (hyperparam_fname, 'w') as hp_f : 
        for combinaison in lbd_eta_mu_values : 
            hp_f.write("%s\n" % combinaison)
    #-----------------------------------

    #-----------------------------------
    # For each fold : 
    # use subsamble to test feature selection with combinaisons of hyperparameters from the grid
    # find opt param
    # use these opt param to select feature using the all training set. = predict causal status of features <- quantify these perf
    # use a ridge regression trained with selected features only to predict quantitativ phenotypes on test set <- quantify these perf
    if SEQ_MODE : 
        for fold_idx in xrange(args.num_folds):
            logging.info ("============= FOLD : %d"%fold_idx)
            run_fold(
                fold_idx,
                args, 
                lbd_eta_values, lbd_eta_mu_values_np, lbd_eta_mu_values, 
                evalf.xp_indices[fold_idx], 
                genotype_fname, network_fname , precision_fname , causal_fname, phenotype_fnames, scores_fnames,
                resu_dir)
            run_predictions(fold_idx, args, resu_dir, data_dir, ef.xp_indices[fold_idx]['trIndices'], ef.xp_indices[fold_idx]['teIndices']) #XXX here ? or in main ?
        # END for fold_idx in range(args.num_folds)
    else :
        cmd = "qsub -cwd -V -N r%df -t 1-%d \
               qsub_run-fold.sh  %d %d %d %d %d %d %s %s %s %s %d" \
               %( 
                  repeat_idx, args.num_folds,
                  args.num_tasks, args.num_features, args.num_samples, args.num_repeats, args.num_folds, args.num_subsamples,
                  args.data_dir, args.resu_dir, args.simu_id, hyperparam_fname, repeat_idx)
        print cmd
        p = subprocess.Popen(shlex.split(cmd))
        # run predictions -> in main
        # END for fold_idx in range(args.num_folds)
        

    if SEQ_MODE : 
        print_analysis_files(args, resu_dir, data_dir,  evalf.xp_indices)
    else : 
        cmd = "qsub -cwd -V -hold_jid r%df ./qsub_print_analysis_files.sh\
            %d %d %d %d %d %d %s %s %s %d" \
            %( repeat_idx, 
               args.num_tasks, args.num_features, args.num_samples, args.num_repeats, args.num_folds, args.num_subsamples,
               args.data_dir, args.resu_dir, args.simu_id, repeat_idx)
        p = subprocess.Popen(shlex.split(cmd))

def run_predictions(fold_idx, args, resu_dir, data_dir, trIndices, teIndices ):
    # can't be in the qsub run fold due to pytable utilisation.

    genotype_fname = '%s/%s.genotypes.txt' % (data_dir, args.simu_id)
    phenotype_fnames = ['%s/%s.phenotype_%d.txt' % \
                        (data_dir, args.simu_id, task_idx) \
                        for task_idx in range(args.num_tasks)]

    #------------------------------------------------------------------
    logging.info( "======== Prediction using opt param")

    # need to know selected_<algo> (list of list :list of sel features, for each task)
    # saved in files '<resu_dir>/<args.simu_id>.<algo>.fold_<fold_idx>.selected_features'
    # one line per task, list of sel feature on each line. 
    selected_st = []
    fname = '%s/%s.sfan.fold_%d.selected_features' % \
        (resu_dir, args.simu_id, fold_idx)
    with open(fname, 'r') as f :
        for line in f : #list of selected feature for a task
            selected_st.append([int(x) for x in line.split()])

    selected_np=[]
    fname = '%s/%s.msfan_np.fold_%d.selected_features' % \
        (resu_dir, args.simu_id, fold_idx)
    with open(fname, 'r') as f :
        for line in f : #list of selected feature for a task
            selected_np.append([int(x) for x in line.split()])

    selected = []
    fname = '%s/%s.msfan.fold_%d.selected_features' % \
        (resu_dir, args.simu_id, fold_idx)
    with open(fname, 'r') as f :
        for line in f : #list of selected feature for a task
            selected.append([int(x) for x in line.split()])



    # For each algorithm, for each task,
    # predict on the test set using a ridge-
    # regression trained with the selected features only.

    for task_idx in range(args.num_tasks):
        logging.info ('task n. %d' %task_idx)
        # Single task
        logging.info("st")
        fname = '%s/%s.sfan.fold_%d.task_%d.predicted' % \
                (resu_dir, args.simu_id, fold_idx, task_idx)
        ef.run_ridge_selected(selected_st[task_idx], genotype_fname,
                              phenotype_fnames[task_idx],
                              trIndices, teIndices, fname)

        # Multitask (no precision)
        logging.info("np")
        fname = '%s/%s.msfan_np.fold_%d.task_%d.predicted' % \
                (resu_dir, args.simu_id, fold_idx, task_idx)
        ef.run_ridge_selected(selected_np[task_idx], genotype_fname,
                              phenotype_fnames[task_idx],
                              trIndices, teIndices, fname)

        # Multitask (precision)
        logging.info("msfan")
        fname = '%s/%s.msfan.fold_%d.task_%d.predicted' % \
                (resu_dir, args.simu_id, fold_idx, task_idx)
        ef.run_ridge_selected(selected[task_idx], genotype_fname,
                              phenotype_fnames[task_idx],
                              trIndices, teIndices, fname)
    #------------------------------------------------------------------

def print_analysis_files(args, resu_dir, data_dir, xp_indices):
    """TODO
    """

    analysis_files = get_analysis_files_names(args.resu_dir, args.simu_id)

    phenotype_fnames = ['%s/%s.phenotype_%d.txt' % \
                        (data_dir, args.simu_id, task_idx) \
                        for task_idx in range(args.num_tasks)]


    # Concatenate ppv and tpr files :
    bash_cmd = "head -c -1 -q"
    #bash_cmd = "cat %s/repeat%d/%s.<algo>.fold*.ppv | tr -d '\n'" % (args.resu_dir, repeat_idx, args.simu_id)
    # problem : 
    # Expanding the * glob is part of the shell, 
    # but by default subprocess does not send your commands via a shell, 
    # so the command is executed (ls, head, or whatever), then a literal * is used as an argument.
    # -> supply shell=True to execute the command through a shell interpreter
    #    (security priblem)
    # -> use the glob module
    #    (quicker (no process startup overhead), more reliable and cross platform)

    #----------------------------------------------------------------------
    # print to file ppv : 
    logging.info( "======== Print ppv file")
    # for each repeat, each fold print ppv in its own file 
    # we want to concatenate these files in one line we print in analysis_files['pvv<algo>']
       
    template_f_names = resu_dir+"/"+args.simu_id+".%s.fold_*.ppv"

    process = subprocess.Popen(shlex.split(bash_cmd)+glob.glob(template_f_names % 'sfan'), stdout=subprocess.PIPE)
    output = process.communicate()[0]
    with open(analysis_files['ppv_st'], 'a') as f:
        f.write("%s\n" %output)
    
    process = subprocess.Popen(shlex.split(bash_cmd)+glob.glob(template_f_names % 'msfan_np'), stdout=subprocess.PIPE)
    output = process.communicate()[0]
    with open(analysis_files['ppv_msfan_np'], 'a') as f:
        f.write("%s\n" %output)

    process = subprocess.Popen(shlex.split(bash_cmd)+glob.glob(template_f_names % 'msfan'), stdout=subprocess.PIPE)
    output = process.communicate()[0]
    with open(analysis_files['ppv_msfan'], 'a') as f:
        f.write("%s\n" %output)

    #----------------------------------------------------------------------
    # print to file tpr : 
    logging.info( "======== Print tpr file")
    # for each repeat, each fold print pr in its own file 
    # we want to concatenate these files in one line we print in analysis_files['tpr<algo>']
    template_f_names = resu_dir+"/"+args.simu_id+".%s.fold_*.tpr"
    
    process = subprocess.Popen(shlex.split(bash_cmd)+glob.glob(template_f_names % 'sfan'), stdout=subprocess.PIPE)
    output = process.communicate()[0]
    with open(analysis_files['tpr_st'], 'a') as f:
        f.write("%s\n" %output)

    process = subprocess.Popen(shlex.split(bash_cmd)+glob.glob(template_f_names % 'msfan_np'), stdout=subprocess.PIPE)
    output = process.communicate()[0]
    with open(analysis_files['tpr_msfan_np'], 'a') as f:
        f.write("%s\n" %output)

    process = subprocess.Popen(shlex.split(bash_cmd)+glob.glob(template_f_names % 'msfan'), stdout=subprocess.PIPE)
    output = process.communicate()[0]
    with open(analysis_files['tpr_msfan'], 'a') as f:
        f.write("%s\n" %output)

                        
    #----------------------------------------------------------------------
    logging.info( "======== Compute RMSE")
    # For each algorithm, and for each task, compute RMSE
    # using :
    #   - an external function 
    #   - the predictions saved in files (fold per fold)
    #   - the true values given by phenotype_fnames[task_idx] and te_indices
    # save to file '%s/%s.<algo>.rmse' % (args.resu_dir, args.simu_id)
    # => rmse_st_fname ; rmse_np_fname ; rmse_fname
    # Files structure : 
    # each line = a repeat
    # on each line there are several RMSE values, one per task

    # Single task
    predicted_fname = resu_dir+'/'+args.simu_id+'.sfan.fold_%d.task_%d.predicted' 
    rmse_list = ef.compute_ridge_selected_RMSE( phenotype_fnames, predicted_fname, 
                                    xp_indices, args.num_tasks)
    with open(analysis_files['rmse_st'], 'a') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in rmse_list]))
    # Multitask (no precision)
    predicted_fname = resu_dir+'/'+args.simu_id+'.msfan_np.fold_%d.task_%d.predicted' 
    rmse_list = ef.compute_ridge_selected_RMSE( phenotype_fnames, predicted_fname, 
                                    xp_indices, args.num_tasks)
    with open(analysis_files['rmse_msfan_np'], 'a') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in rmse_list]))
    # Multitask (precision)
    predicted_fname = resu_dir+'/'+args.simu_id+'.msfan.fold_%d.task_%d.predicted' 
    rmse_list = ef.compute_ridge_selected_RMSE( phenotype_fnames, predicted_fname, 
                                    xp_indices, args.num_tasks)             
    with open(analysis_files['rmse_msfan'], 'a') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in rmse_list]))
    #----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    logging.info( "======== Compute CI ")
    # For each algorithm, and for each task, compute consistency index
    # between the features selected for each fold.
    # Use an external function using ef.consistency_index_k()
    # use the selected features saved to files and the true causal features
    # save to file '%s/%s.<algo>.consistency' % (args.resu_dir, args.simu_id)
    # File structure : 
    # each line = a repeat
    # on each line there are several ci values, one per task

    # Single task
    selection_fname = resu_dir+'/'+args.simu_id+'.sfan.fold_%d.selected_features'
    ci_list = ef.consistency_index_task(selection_fname, args.num_folds, args.num_tasks, args.num_features)
    with open(analysis_files['ci_st'], 'a') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in ci_list]))
    # Multitask (no precision)
    selection_fname = resu_dir+'/'+args.simu_id+'.msfan_np.fold_%d.selected_features'
    ci_list = ef.consistency_index_task(selection_fname, args.num_folds, args.num_tasks, args.num_features)
    with open(analysis_files['ci_msfan_np'], 'a') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in ci_list]))
    # Multitask (precision)
    selection_fname = resu_dir+'/'+args.simu_id+'.msfan.fold_%d.selected_features'
    ci_list = ef.consistency_index_task(selection_fname, args.num_folds, args.num_tasks, args.num_features)
    with open(analysis_files['ci_msfan'], 'a') as f:
        f.write('%s \n' % ' '.join(['%.2f ' % x for x in ci_list]))

    #-----------------------------------------------------------------------
    logging.info( "analysis_files outputed")




def print_plot_files(f_name, means_to_plot, std_to_plot):
    """ Print some means and std in a plot module understandable file format.

    Parameters
    ----------
    f_name : filename
        Path to the output file.
    means_to_plot : list of means 
        (one per task)
    std_to_plot : list of std
        (one per task)

    /!\ : 
    len(means_to_plot) and len(std_to_plot) should be equal.

    Side effects
    ------------
    Print a file named <f_name>. 
    Its format is the following : 
    Space-separated list of means (one per task) | Space-separated list of std (one per task)
    (one line per algo : st, np and msfan)
    """
    algos = ['st', 'np', 'msfan']
    with open(f_name, 'w') as f: 
        for algo in algos : 
            f.write(
                ' '.join (str(item) for item in means_to_plot[algo])+
                '|'+
                ' '.join (str(item) for item in std_to_plot[algo])+'\n'
            )

def print_save_res_measures(means, std, fname):
    """ Print out & save with LaTeX format a measure table.
    
    Parameters
    ----------
    means : list of means
        (one per measure)
    std : list of std
        (one per measure)
    fname : filename
        Path to the output file.

    /!\ : 
    means and std should have the same length
    
    Side effects
    ------------
    Create a file named <fname> holding the table in LaTeX format.

    """
    #------------------
    # Print out measures tables
    # and save them (with LaTeX table format) in plain text file

    header_print = (
            "-----------------------------------------------------------------------\n"
            "{:^80}\n"
            "       +--------------------------------------------------------------+\n"
            "       |                              algo                            |\n"
            "  task |         sfan       |          np        |          msfan     |\n"
            "=======+====================+====================+====================+\n"
    )
    header_save = "task & sfan &np &msfan \\\\\\hline\n"
    algos = ['st', 'np', 'msfan']
    for measure in means : 
        to_print = ''
        to_save = ''
        for task_idx in xrange (len (means[measure]['st'] ) ) : 
            to_print += '{:^7d}|'.format(task_idx) 
            to_save += '{:^7d}&'.format(task_idx) 
            for algo in algos : 
                to_print += '{:9.3f} ±{:9.3f}|'.format(means[measure][algo][task_idx], std[measure][algo][task_idx]) 
                to_save += '{:9.3f} ±{:9.3f}&'.format(means[measure][algo][task_idx], std[measure][algo][task_idx]) 
            to_save = to_save[:-1] #don't want the last '&'
            to_save += "\\\\\n" # two step and not to_save[-1] = "\\\\\n" because python strings are immuable
            to_print+="\n"
        
        print header_print.format(measure) + to_print

        with open(fname%measure, 'w') as f : 
            f.write(header_save+to_save)



def handle_measures_results(analysis_files, args): 
    """TODO
    """
    # For each measure compute average/mean +- standard deviation per task for 
    means, std = extract_res_means_and_std(analysis_files, args)
    

    fname = args.resu_dir+'/'+args.simu_id+'.results_%s'
    print_save_res_measures(means, std, fname)

    # Plots : 

    for measure in means : 
        f_name = "%s/%s.%s_plot.values" %(args.resu_dir, args.simu_id, measure)
        print_plot_files(f_name, means[measure], std[measure])
        plot.bar_plot(measure, f_name) 


def extract_res_means_and_std(analysis_files, args):
    """
    TODO
    """
    # For each measure compute average/mean +- standard deviation per task for each algo
    means = {}
    std = {}

    # for : 

    # RMSE : 
    means["rmse"], std["rmse"] = ef.extract_res_from_files(
        [   analysis_files['rmse_st'],
            analysis_files['rmse_msfan_np'],
            analysis_files['rmse_msfan']
        ],
        args.num_tasks,
        args.num_repeats
    ) 
    # consistency index : 
    means["ci"], std["ci"]= ef.extract_res_from_files(
        [   analysis_files['ci_st'],
            analysis_files['ci_msfan_np'],
            analysis_files['ci_msfan']
        ],
        args.num_tasks,
        args.num_repeats
    )
    # PPVs : 
    means["ppv"], std["ppv"] = ef.extract_res_from_files(
        [   analysis_files['ppv_st'],
            analysis_files['ppv_msfan_np'],
            analysis_files['ppv_msfan']
        ],
        args.num_tasks,
        args.num_repeats,
        args.num_folds # needed to handle particular file organisation fer foldt then per task
    )
    # sensitivities : 
    means["tpr"],std["tpr"] = ef.extract_res_from_files(
        [   analysis_files['tpr_st'],
            analysis_files['tpr_msfan_np'],
            analysis_files['tpr_msfan']
        ],
        args.num_tasks,
        args.num_repeats,
        args.num_folds # needed to handle particular file organisation fer foldt then per task
    )
    return means, std

def handle_measures_results(analysis_files, args): 
    """TODO
    """
    # For each measure compute average/mean +- standard deviation per task for each algo
    means, std = extract_res_means_and_std(analysis_files, args)
    

    fname = args.resu_dir+'/'+args.simu_id+'.results_%s'
    print_save_res_measures(means, std, fname)

    # Plots : 

    for measure in means : 
        f_name = "%s/%s.%s_plot.values" %(args.resu_dir, args.simu_id, measure)
        print_plot_files(f_name, means[measure], std[measure])
        if SEQ_MODE : #XXX cluster  no display name and no $DISPLAY environment variable
            plot.bar_plot(measure, f_name) 


def main():
    """ Sequentially run validation experiments on synthetic data.

    Arguments
    ---------
    args.num_tasks: int
        Number of tasks.
    args.num_features: int
        Number of features.
    args.num_samples: int
        Number of samples
    args.num_repeats: int
        Number of repeats.
    args.num_folds: int
        Number of cross-validation folds.
    args.num_subsamples: int
        Number of subsamples (for stability evaluation).
    args.data_dir: filename
        Path of the directory in which to save the simulated data.
    args.resu_dir: filename
        Path of the directory in which to save the simulated data.
    args.simu_id: string
        Name of the simulation, to be used to name files within args.root_dir.
    args.verbose: boolean
        If true, turn on detailed information logging.

    Generated files
    ---------------
    1. Simulated data
    For each repeat, under <data_dir>/repeat_<repeat_id>:
        <simu_id>.readme:
            README file describing the simulation paramters
        <simu_id>.task_similarities.txt:
            Matrix of precision between tasks
        <simu_id>.causal_features:
            Lists of causal features.
            One list per task. Indices start at 0.
        <simu_id>.causal_weights:
            Lists of the weights given to the causal features.
            One list per task. Same order as in <data_dir>/<simu_id>.causal_features.
        <simu_id>.genotypes.txt:
            num_features x num_samples matrix of {0, 1, 2} (representing SNPs).
        <simu_id>.network.dimacs:
            Network over the features.
        For task_id in 0, ..., args.num_tasks:
            <simu_id>.phenotype_<task_id>.txt:
                Phenotype vector (of size args.num_samples) for task <task_id>.
            <simu_id>.scores_<task_id>.txt
                Node weights (of size args.num_features) for task <task_id>.
                Computed as Pearson correlation.

        For each fold_idx:
            <simu_id>.<fold_idx>.trIndices
                Space-separated list of training indices.
            <simu_id>.<fold_idx>.teIndices
                Space-separated list of test indices.
            For each subsample_idx:
                <data_dir>/<simu_id>.<fold_idx>.<ss_idx>.ssIndices
                    Space-separated lists of subsample indices,
                    one line per list / subsample.

    2. Results
    For each repeat, under <resu_dir>/repeat_<repeat_idx>:
        For each fold_idx:
            For each algo in ('sfan', 'msfan_np', 'msfan'):
                <simu_id>.<algo>.fold_<fold_idx>.parameters
                     Optimal parameters.

                <simu_id>.<algo>.fold_<fold_idx>.selected_features
                    List of list of selected features, one per task.
                    Each line corresponds to a task and contains a 
                    space-separated list of indices, STARTING AT 0.

                <simu_id>.<algo>.fold_<fold_idx>.task_<task_idx>.predicted
                    Predictions, on the test set, of a ridge regression
                    trained only on the selected features.
    
    Under <resu_dir>:
        For each algo in ('sfan', 'msfan_np', 'msfan'):
            <simu_id>.<algo>.rmse:
                Space-separated lists of RMSEs (one per task)
                each line corresponds to one repeat. 
            <simu_id>.<algo>.consistency:
                Space-separated lists of consistency index (one per task)
                each line corresponds to one repeat. 
            <simu_id>.<algo>.ppv
                Space-separated lists of PPVs per task, per fold,
                each line corresponds to one repeat. 
            <simu_id>.<algo>.sensitivity
                Space-separated lists of sensitivities per task, per fold,
                each line corresponds to one repeat. 
            <simu_id>.results
                Average/standard deviation values for: consistency index, RMSE,
                PPV and sensitivity.
            TODO: Plot files.
    
    For file format specifications, see README.md


    
    Example
    -------
    $ python synthetic_data_experiments.py -k 3 -m 200 -n 100 -r 10 -f 10 -s 10 \
             ../data/simu_synth_01 ../results/simu_synth_01 simu_01 --verbose


    ### Sur le cluster : XXX

    $ python synthetic_data_experiments.py -k 3 -m 50 -n 35 -r 5 -f 3 -s 2 \
             /share/data40T/athenais/data/simu_synth_01 /share/data40T/athenais/results/simu_synth_01 simu_01 --verbose
    """

    # Arguments handling : 
    args = get_integrous_arguments_values()



    #-------------------------------------------------------------------------
    # Insure working repositories exist : 

    # Create simulated data repository if it does not exist
    create_dir_if_not_exists(args.data_dir)

    # (Delete and re)create results repository (if it exists)
    if os.path.isdir(args.resu_dir): 
        logging.info("Deleting %s\n" % args.resu_dir)
        try:
            shutil.rmtree(args.resu_dir)
        except OSError:
            raise #???XXX
    # More pythonic : ??? 
    #import errno
    #try:
    #    os.remove(args.resu_dir)
    #except OSError as e :
    #    if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
    #       raise
    create_dir_if_not_exists(args.resu_dir)
    #-------------------------------------------------------------------------



    #-------------------------------------------------------------------------
    # Create analysis file names :
    analysis_files = get_analysis_files_names(args.resu_dir, args.simu_id)
    #-------------------------------------------------------------------------



    #-------------------------------------------------------------------------
    for repeat_idx in xrange(args.num_repeats):
            run_repeat(repeat_idx, args, analysis_files)

    if not SEQ_MODE : 
        cmd = "python run_predictions.py\
            %d %d %d %d %d %d %s %s %s" \
            %( args.num_tasks, args.num_features, args.num_samples, args.num_repeats, args.num_folds, args.num_subsamples,
            args.data_dir, args.resu_dir, args.simu_id)
        with open('launcher_handle-measures-results.sh', 'a') as f : 
            f.write(cmd)
    """
    if SEQ_MODE : 
        for repeat_idx in xrange(args.num_repeats):
            run_repeat(repeat_idx, args, analysis_files)
    else : 
        # run a job array of which each job is a repeat. 
        # -N    : job array name is 'repeats' 
        # -cwd  : execute the job in the current work directory 
        #         make scritp knowing each other. 
        # -j y  : join stout and stderr
        # -o <path> : Place the joined output in another location other than the working directory XXX
        # -V : pass all environement variables
        # TODO : use -v VAR1="hello",VAR2="Sreedhar",VAR3="How are you?" to pass variables ? 
        cmd = "qsub -cwd -V -t 1-%d -N repeats\
               qsub_run-repeat.sh  %d %d %d %d %d %d %s %s %s" \
               %( args.num_repeats,
                  args.num_tasks, args.num_features, args.num_samples, args.num_repeats,args.num_folds, args.num_subsamples,
                  args.data_dir, args.resu_dir, args.simu_id)
        print cmd
        p = subprocess.Popen(shlex.split(cmd))
        # Even though shlex.split strips off the single quotes, it appears to be interpreted fine when using subprocess.

    # END for repeat_idx in range(args.num_repeats)
    #-------------------------------------------------------------------------
   
    if SEQ_MODE : 
        #-------------------------------------------------------------------------
        # Handle measures results : 
        handle_measures_results(analysis_files, args)
        #------------------
    else : 
        # run this job when repeats are finished. 
        # -hold_jid repeats : hold job starting until the job named repeat is complete. 
        # -m e -M vaginay.athenais@gmail.com : send me an email when the job is finished. 
        # -V : pass all environement variables
        # TODO : use -v VAR1="hello",VAR2="Sreedhar",VAR3="How are you?" to pass variables ? 
        cmd = "qsub -cwd -V -hold_jid repeats ./qsub_handle-measures-results.sh\
            %d %d %d %d %d %d %s %s %s" \
            %( args.num_tasks, args.num_features, args.num_samples, args.num_repeats, args.num_folds, args.num_subsamples,
            args.data_dir, args.resu_dir, args.simu_id)
        p = subprocess.Popen(shlex.split(cmd))
    """

    cmd = "qsub -cwd -V ./qsub_handle-measures-results.sh\
            %d %d %d %d %d %d %s %s %s" \
            %( args.num_tasks, args.num_features, args.num_samples, args.num_repeats, args.num_folds, args.num_subsamples,
            args.data_dir, args.resu_dir, args.simu_id)
    with open('launcher_handle-measures-results.sh', 'a') as f : 
        f.write(cmd)

if __name__ == "__main__":
    
    main()
    
    print("THE END")
          
