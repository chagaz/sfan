"""synthetic_data_experiments.py -- Run validation experiments on synthetic data

In this version all experiments are run sequentially.
"""
# TODO: Deal with the fact that we're doing all of the evaluation
# for each algorithm; maybe this should be a separate function.

# TODO: Separate work in more functions.
    
# Also deal with the fact that we have different metrics for each task,
# resultint in a different structure for the results files.
 



import argparse
import logging
import os
import subprocess
import sys

import evaluation_framework
import generate_data
import multitask_sfan


def main():
    """
    Sequentially run validation experiments on synthetic data.

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
            Matrix of correlation between tasks
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
            <simu_id>.<fold_idx>.selected_features
                Space-separated list of indices of selected_features,
                STARTING AT 0.
            <simu_id>.<fold_idx>.predicted
                Predicted phenotypes for the test set, one value per line.
            <simu_id>.<fold_idx>.parameters
                Optimal parameters for each of single-task,
                multi-task with no correlation and multi-task with correlation.

    Under <resu_dir>:
        <simu_id>.rmse:
            List of final RMSEs (one per repeat),
            one value per line.
        <simu_id>.consistency:
            List of final Consistency Indices (one per repeat),
            one value per line.
        <simu_id>.ppv
            Space-separated lists of PPVs (one value per fold),
            each line corresponds to one repeat. 
        <simu_id>.sensitivity
            Space-separated lists of sensitivities (one value per fold),
            each line corresponds to one repeat. 
        <simu_id>.results
            Average/standard deviation values for: consistency index, RMSE,
            PPV and sensitivity.
    
    For file format specifications, see README.md


    
    Example
    -------
    $ python synthetic_data_experiments -k 3 -m 200 -n 100 -r 10 -f 10 -s 10 \
             ../data/simu_synth_01 ../results/simu_synth_01 simu_01 --verbose
    """
    # Get arguments values
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
    args = parser.parse_args()

    # Check arguments integrity
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

    # Create simulated data repository if it does not exist
    if not os.path.isdir(args.data_dir):
        logging.info("Creating %s\n" % args.data_dir)
        try: 
            os.makedirs(args.data_dir)
        except OSError:
            if not os.path.isdir(args.data_dir):
                raise

    # Create results repository if it does not exist
    if not os.path.isdir(args.resu_dir):
        logging.info("Creating %s\n" % args.resu_dir)
        try: 
            os.makedirs(args.resu_dir)
        except OSError:
            if not os.path.isdir(args.resu_dir):
                raise


    for repeat_idx in range(args.num_repeats):
        # Instantiate data generator
        data_dir = '%s/repeat_%d' % (args.data_dir, repeat_idx)

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
        network_fname = '%s/%s.network.dimacs' % (data_dir, args.simu_id)
        correlation_fname = '%s/%s.task_similarities.txt' % (data_dir,
                                                             args.simu_id)
        phenotype_fnames = ['%s/%s.phenotype_%d.txt' % \
                            (data_dir, args.simu_id, task_idx) \
                            for task_idx in range(num_tasks)]

        

        # Instantiate evaluation framework
        ef = evaluation_framework.Framework(args.num_samples, args.num_folds,
                                            args.num_subsamples)

        # Compute cross-validation folds and subsample indices
        ef.compute_indices()

        # Save cross-validation folds and subsample indices to file
        ef.save_indices(data_dir, args.simu_id)


        # TODO: Create <resu_dir>/repeat_<repeat_id> if it does not exist
        

        
        # TODO (CA): define the ranges of values for lbd, eta, mu
        # Note: acceptable values of lbd, eta and mu are linked,
        # we won't define all three ranges as constants.

        # Note: MU_VALUES must contain 0, as this corresponds to
        # solving all tasks separately (single-task)
        lbd_values = [] # <-- TODO 
        eta_values = [] # <-- TODO 
        mu_values = [0] # <-- TODO 
        
        for fold_idx in range(args.num_folds):
            # Inititalize dictionary to store selected features
            # sf is a list of lists of selected features
            # (one per subsample iteration)
            # sf_dict is a nested dictionary, indexed by
            #   - value of mu
            #   - value of lbd
            #   - value of eta
            #   - value of task_idx
            # i.e. you get a specific sf from sf_dict by querying
            # sf_dict[str(mu)][str(lbd)][str(eta)][task_idx]
            sf_dict = {}          # using correlation matrix
            sf_nc_dict = {}       # not using correlation matrix

            # Use string representations of lbd, eta, mu
            lbd_s = str(lbd)
            eta_s = str(lbd)
            mu_s = str(lbd)
            
            for lbd in range(lbd_values):
                sf_dict[lbd_s] = {}          
                sf_nc_dict[lbd_s] = {}
                for eta in range(eta_values):
                    sf_dict[lbd_s][eta_s] = {}          
                    sf_nc_dict[lbd_s][eta_s] = {}
                    for mu in range(mu_values):
                        sf_dict[lbd_s][eta_s][mu_s] = {}          
                        sf_nc_dict[lbd_s][eta_s][mu_s] = {}
                        for task_idx in range(args.num_tasks):
                            sf_dict[lbd_s][eta_s][mu_s][task_idx] = []
                            sf_nc_dict[lbd_s][eta_s][mu_s][task_idx] = [] 
                
            for ss_idx in range(args.num_subsamples):
                # Get samples
                sample_indices = ef.xp_indices[fold_idx]['bsIndices'][ss_idx]

                # TODO:
                # Generate sample-specific network scores
                # Save to files
                weights_fname = []

                
                for lbd in range(lbd_values):
                    lbd_s = str(lbd)
                    for eta in range(eta_values):
                        eta_s = str(lbd)
                        for mu in range(mu_values):
                            mu_s = str(lbd)

                            sf_nc_d = sf_nc_dict[mu_s][lbd_s][eta_s]
                            sf_d = sf_dict[mu_s][lbd_s][eta_s]
                            
                            # Run multi-task sfans *without* correlation matrix
                            # Ideally, I'd do the following:
                            # sfan_solver = Sfan(args.num_tasks, network_fname,
                            #                    weights_fname,
                            #                    lbd, eta, mu,
                            #                    correlation_fname,
                            #                    output_f='/tmp/test')
                            # tt = sfan_solver.create_dimacs()
                            # sfan_solver.run_maxflow()

                            # But because cython output to screen is NOT
                            # caught by sys.stdout,
                            # we need to run this externally...
                            p = subprocess.Popen(['python',
                                                  'multitask_sfan.py',
                                                  '--num_tasks',
                                                  str(args.num_tasks),
                                                  '--networks', network_fname,
                                                  '--node_weights',
                                                  weights_fnames,
                                                  '-l', lbd_s, '-e', eta_s,
                                                  '-m', mu_s,
                                                  '--output', '/tmp/test'], 
                                                  stdout=subprocess.PIPE)
                            p_out = p.communicate()[0].split("\n")[2:]

                            # Process the output to get lists of selected
                            # features
                            sel_ = [[int(x) for x in line.split()] \
                                         for line in p_out]

                            # Store selected features in the dictionary
                            for task_idx, sel_list in enumerate(sel_):
                                sf_nc_d[task_idx].append(sel_list)

                            # run multi-task sfans *with* correlation matrix
                            p = subprocess.Popen(['python',
                                                  'multitask_sfan.py',
                                                  '--num_tasks',
                                                  str(args.num_tasks),
                                                  '--networks', network_fname,
                                                  '--node_weights',
                                                  weights_fnames,
                                                  '--correlation_matrix',
                                                  correlation_fname,
                                                  '-l', lbd_s, '-e', eta_s,
                                                  '-m', mu_s,
                                                  '--output', '/tmp/test'], 
                                                  stdout=subprocess.PIPE)
                            p_out = p.communicate()[0].split("\n")[2:]

                            # Process the output to get lists of selected
                            # features
                            sel_ = [[int(x) for x in line.split()] \
                                         for line in p_out]

                            # Store selected features in the dictionary
                            for task_idx, sel_list in enumerate(sel_):
                                sf_d[task_idx].append(sel_list)
            # END for ss_idx in range(args.num_subsamples)
                            
            # Compute consistencies for single-task,
            # multi-task with no correlation matrix,
            # multi-task with correlation matrices.

            # Get optimal parameter values for each algo.
            # single task: 
            lbd_opt_st = 0
            eta_opt_st = 0

            # multitask, no correlation: 
            lbd_opt_nc = 0
            eta_opt_nc = 0
            mu_opt_nc = 0
            
            # multitask, with correlation: 
            lbd_opt = 0
            eta_opt = 0
            mu_opt = 0
            
            for lbd in range(lbd_values):
                for eta in range(eta_values):

                    # Single-task: mu = 0
                    # TODO: Compute average consistency index across tasks:
                    for task_idx in range(args.num_tasks):
                        sel_list = sf_dict[str(0.)][lbd_s][eta_s][task_idx]
                        cidx = evaluation_framework.consistency_index_k(sel_list)
                        

                    # TODO: Update lbd_opt_st, eta_opt_st accordingly

                    # Multi-task
                    # TODO Exclude mu=0 from mu_values when looping
                    for mu in range(mu_values):
                        # Multitask with no correlation
                        sel_list_d = sf_nc_dict[mu_s]
                        for task_idx in range(args.num_tasks):
                            # TODO: Compute average consistency index across tasks:
                            sel_list = sel_list_d[lbd_s][eta_s][task_idx]
                            cidx = evaluation_framework.consistency_index_k(sel_list)
                        # TODO: Update lbd_opt_nc, eta_opt_nc, mu_opt_nc
                        # accordingly

                        # Multitask with correlation
                        sel_list_d = sf_dict[mu_s]
                        for task_idx in range(args.num_tasks):
                            # TODO: Compute average consistency index across tasks:
                            sel_list = sel_list_d[lbd_s][eta_s][task_idx]
                            cidx = evaluation_framework.consistency_index_k(sel_list)
                        # TODO: Update lbd_opt, eta_opt, mu_opt
                        # accordingly

            # TODO: For each algorithm, save optimal parameters to file
            fname = '%s/repeat_%d/%s/%d.parameters' % (args.resu_dir,
                                                      repeat_idx, simu_id,
                                                      fold_idx)
                            
            # TODO: For each algorithm, run algorithms again to select features,
            # using the whole training set

    
            # TODO: For each algorithm, save selected features to file
            fname = '%s/repeat_%d/%s/%d.selected_features' % (args.resu_dir,
                                                              repeat_idx, simu_id,
                                                              fold_idx)
                            
            # TODO: For each algorithm, train a ridge-regression
            # on the selected features / training set
            # Use sklearn.linear_models


            # TODO: For each algorithm, predict on test set


            # TODO: For each algorithm, save predictions to file
            fname = '%s/repeat_%d/%s/%d.predicted' % (args.resu_dir,
                                                      repeat_idx, simu_id,
                                                      fold_idx)
            

            # TODO: For each algorithm, and for each task,
            # compute PPV and sensitivity,
            # comparing to <root_dir>/<simu_id>.causal_features

            # Use sklearn.metrics

        
            # TODO: Save PPV to file (use 'a' mode)
            fname = '%s/%s.ppv' % (args.resu_dir, simu_id)
            # TODO: Save sensitivity to file (use 'a' mode):
            fname = '%s/%s.sensitivity' % (args.resu_dir, simu_id)
            
                            
        # END for fold_idx in range(args.num_folds)

                            
    # TODO: Compute RMSE (use sklearn.metrics)


    # TODO: Save RMSE to file
    fname = '%s/%s.rmse' % (args.resu_dir, simu_id)
            

    # TODO: Compute consistency index


    # TODO: Save consistency index to file
    fname = '%s/%s.consistency' % (args.resu_dir, simu_id)


    # END for repeat_idx in range(args.num_repeats)



                            
    # TODO: Compute average/mean for RMSE, PPVs, sensitivities, CI.

    
    # TODO: Print out and save in plain text file
    fname = '%s/%s.results' % (args.resu_dir, simu_id)


    # TODO: Display (use matplotlib)
                            


        