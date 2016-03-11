"""synthetic_data_experiments.py -- Run validation experiments on synthetic data

In this version all experiments are run sequentially.
"""

# Importing local libraries first,
# because otherwise Error in `python': free(): invalid pointer
import multitask_sfan
import evaluation_framework as ef
import generate_data


import argparse
import logging
import os
import numpy as np
import scipy.stats as st
import subprocess
import sys
import tables as tb
import tempfile


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
                List of final RMSEs (one per repeat),
                one value per line.
            <simu_id>.<algo>.consistency:
                List of final Consistency Indices (one per repeat),
                one value per line.
            <simu_id>.<algo>.ppv
                Space-separated lists of PPVs (one value per task and per fold),
                each line corresponds to one repeat. 
            <simu_id>.<algo>.sensitivity
                Space-separated lists of sensitivities (one value per task and per fold),
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

    # Create files to hold PPV and sensitivity values
    ppv_st_fname = '%s/%s.sfan.ppv' % (args.resu_dir, args.simu_id)
    tpr_st_fname = '%s/%s.sfan.sensitivity' % (args.resu_dir, args.simu_id)

    ppv_np_fname = '%s/%s.msfan_np.ppv' % (args.resu_dir, args.simu_id)
    tpr_np_fname = '%s/%s.msfan_np.sensitivity' % (args.resu_dir, args.simu_id)

    ppv_fname = '%s/%s.msfan.ppv' % (args.resu_dir, args.simu_id)
    tpr_fname = '%s/%s.msfan.sensitivity' % (args.resu_dir, args.simu_id)

    # TODO: Create files to hold RMSE and consistency values

    for repeat_idx in range(args.num_repeats):
        print "==================================================== REPETITION :"+`repeat_idx`
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

        

        # Instantiate evaluation framework
        evalf = ef.Framework(args.num_samples, args.num_folds,
                             args.num_subsamples)

        # Compute cross-validation folds and subsample indices
        evalf.compute_indices()

        # Save cross-validation folds and subsample indices to file
        evalf.save_indices(data_dir, args.simu_id)


        # Create <resu_dir>/repeat_<repeat_id> if it does not exist
        resu_dir = "%s/repeat_%d" % (args.resu_dir, repeat_idx)
        if not os.path.isdir(resu_dir):
            logging.info("Creating %s\n" % resu_dir)
            try:
                os.makedirs(resu_dir)
            except OSError:
                if not os.path.isdir(resu_dir):
                    raise

        #-----------------------------------
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
                                    tmp_scores_f_list, 0, 0, 0,
                                    precision_matrix_f=precision_fname)
        lbd_eta_mu_values = sfan_.compute_hyperparameters_range(num_values=5)
        lbd_eta_values = [" ".join(plist.split()[:-2]) \
                          for plist in lbd_eta_mu_values]

        # Delete temporary files from tmp_scores_f_list
        for fname in tmp_scores_f_list:
            os.remove(fname)
        #-----------------------------------

        for fold_idx in range(args.num_folds):
            print "==================================================== REPETITION :"+`repeat_idx`+"FOLD :"+`fold_idx`
            # Inititalize dictionary to store selected features
            # sf is a list of lists of selected features
            # (one per subsample iteration)
            # sf_dict is a nested dictionary, indexed by
            #   - value of the parameters
            #   - value of task_idx
            # i.e. you get a specific sf from sf_dict by querying
            # sf_dict[params][task_idx]
            sf_st_dict = {}       # single task
            sf_np_dict = {}       # not using precision matrix
            sf_dict = {}          # using precision matrix
            for params in lbd_eta_values:
                sf_st_dict[params] = {}
                for task_idx in range(args.num_tasks):
                    sf_st_dict[params][task_idx] = []
            for params in lbd_eta_mu_values:
                sf_np_dict[params] = {}
                sf_dict[params] = {}
                for task_idx in range(args.num_tasks):
                    sf_np_dict[params][task_idx] = []
                    sf_dict[params][task_idx] = []

            for ss_idx in range(args.num_subsamples):
                print"==================================================== REPETITION :"+`repeat_idx`+"FOLD :"+`fold_idx`+"SS :"+`ss_idx`
                # Get samples
                sample_indices = evalf.xp_indices[fold_idx]['ssIndices'][ss_idx]

                # Generate sample-specific network scores from phenotypes and genotypes
                tmp_weights_f_list = [] # to hold temp files storing these scores
                with tb.open_file(genotype_fname, 'r') as h5f:
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
                    print"\n\n\n==================================================== XXX lbd_eta_values"
                    # Select features with single-task sfan
                    sel_ = ef.run_sfan(args.num_tasks, network_fname,
                                       tmp_weights_f_list, params)
                    if not sel_ : import pdb; pdb.set_trace() 
                    print "SEL_ sf = \n", sel_
                    # Store selected features in the dictionary
                    for task_idx, sel_list in enumerate(sel_):
                        print"==================================================== ___sf_st_dict"
                        sf_st_dict[params][task_idx].append(sel_list)
                    print "\n\nsf_st_dict = \n", sf_st_dict

                for params in lbd_eta_mu_values:
                    print"\n\n\n==================================================== XXX lbd_eta_mu_values"
                    # Select features with multi-task (no precision) sfan
                    sel_ = ef.run_msfan_nocorr(args.num_tasks, network_fname,
                                               tmp_weights_f_list, params)
                    if not sel_ : import pdb; pdb.set_trace()
                    print "SEL_ no corr = \n", sel_
                    # Store selected features in the dictionary
                    for task_idx, sel_list in enumerate(sel_):
                        print"==================================================== ___sf_np_dict"
                        sf_np_dict[params][task_idx].append(sel_list)
                    # Select features with multi-task sfan
                    sel_ = ef.run_msfan(args.num_tasks, network_fname,
                                        tmp_weights_f_list, precision_fname,
                                        params)  
                    print "SEL_ multi_corr = \n", sel_
                    if not sel_ : import pdb; pdb.set_trace()                                       
                    # Store selected features in the dictionary
                    for task_idx, sel_list in enumerate(sel_):
                        print"==================================================== ___ sf_dict"

                        sf_dict[params][task_idx].append(sel_list)
                        
                # Delete the temporary files stored in tmp_weights_f_list
                for fname in tmp_weights_f_list:
                    os.remove(fname)
            # END for ss_idx in range(args.num_subsamples)
                            
            # Get optimal parameter values for each algo.
            # ??? some lists are empty, is it normal ??? 
            opt_params_st = ef.get_optimal_parameters_from_dict(sf_st_dict, args.num_features)
            opt_params_np = ef.get_optimal_parameters_from_dict(sf_np_dict, args.num_features)
            opt_params = ef.get_optimal_parameters_from_dict(sf_dict, args.num_features)

            # For each algorithm, save optimal parameters to file
            print"==================================================== REPETITION :"+`repeat_idx`+"FOLD :"+`fold_idx`+"SS :"+`ss_idx`+"OPT PARAM ALGO SIMPLE"
            # Single task
            fname = '%s/%s.sfan.fold_%d.parameters' % (resu_dir, args.simu_id, fold_idx)
            with open(fname, 'w') as f:
                f.write(opt_params_st)
                f.close()
            print"==================================================== REPETITION :"+`repeat_idx`+"FOLD :"+`fold_idx`+"SS :"+`ss_idx`+"OPT PARAM ALGO MULTISAN"
            # Multitask (no precision)
            fname = '%s/%s.msfan_np.fold_%d.parameters' % (resu_dir, args.simu_id, fold_idx)
            with open(fname, 'w') as f:
                f.write(opt_params_np)
                f.close()
            print"==================================================== REPETITION :"+`repeat_idx`+"FOLD :"+`fold_idx`+"SS :"+`ss_idx`+"OPT PARAM ALGO MULTIAVEC"
            # Multitask (precision)
            fname = '%s/%s.msfan.fold_%d.parameters' % (resu_dir, args.simu_id, fold_idx)
            with open(fname, 'w') as f:
                f.write(opt_params)
                f.close()

            #------------------------------------------------------------------
            # For each algorithm, run algorithms again to select features,
            # (got a list of list : list of selected features for each task)
            # using the whole training set (i.e. scores_fnames)
            # and optimal parameters.

            # DEBUG ???: give arbitrary optimal parameters as there are no returned due to those used at the begining of the experiment 
            opt_params_st = '-l 2 -e 3 -m 1'
            opt_params_np = '-l 2 -e 3 -m 1'
            opt_params = '-l 2 -e 3 -m 1'

            print"==================================================== REPETITION :"+`repeat_idx`+"FOLD :"+`fold_idx`+"SS :"+`ss_idx`+"RUN ALGO SIMPLE"
            selected_st = ef.run_sfan(args.num_tasks, network_fname,
                                       scores_fnames, opt_params_st)
            print"==================================================== REPETITION :"+`repeat_idx`+"FOLD :"+`fold_idx`+"SS :"+`ss_idx`+"RUN ALGO MULTISAN"
            selected_np = ef.run_msfan_nocorr(args.num_tasks, network_fname,
                                               scores_fnames, opt_params_np)
            print"==================================================== REPETITION :"+`repeat_idx`+"FOLD :"+`fold_idx`+"SS :"+`ss_idx`+"RUN PARAM ALGO MULTIAVEC"
            selected = ef.run_msfan(args.num_tasks, network_fname,
                                        scores_fnames, precision_fname,
                                        opt_params)

            # For each algorithm, save selected features to file
            
            # Single task
            fname = '%s/%s.sfan.fold_%d.selected_features' % \
                    (resu_dir, args.simu_id, fold_idx)
            with open(fname, 'w') as f:
                for selected_features_list in selected_st:
                    f.write("%s\n" % ' '.join(str(x) for x in selected_features_list))
                                            # selected_features_list is a list of int 
                                            # that have to be cast as string so we can join them
                f.close()
            # Multitask (no precision)
            fname = '%s/%s.msfan_np.fold_%d.selected_features' % \
                    (resu_dir, args.simu_id, fold_idx)
            with open(fname, 'w') as f:
                for selected_features_list in selected_np:
                    f.write("%s\n" % ' '.join(str(x) for x in selected_features_list))
                f.close()
            # Multitask (precision)
            fname = '%s/%s.msfan.fold_%d.selected_features' % \
                    (resu_dir, args.simu_id, fold_idx)
            with open(fname, 'w') as f:
                for selected_features_list in selected:
                    f.write("%s\n" % ' '.join(str(x) for x in selected_features_list))
                f.close()
            #------------------------------------------------------------------
            

            #-----------------------------------------------------------
            # For each algorithm, and for each task, compute PPV
            # and sensitivity, and save to ppv_fname, tpr_fname
            
            # Single task
            ppv_list, tpr_list = ef.compute_ppv_sensitivity(causal_fname,
                                                            selected_st,
                                                            args.num_features)
            with open(ppv_st_fname, 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in ppv_list]))
                f.close()

            with open(tpr_st_fname, 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in tpr_list]))
                f.close()


            # Multitask (no precision)
            ppv_list, tpr_list = ef.compute_ppv_sensitivity(causal_fname,
                                                            selected_np,
                                                            args.num_features)
            with open(ppv_np_fname, 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in ppv_list]))
                f.close()

            with open(tpr_np_fname, 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in tpr_list]))
                f.close()

            # Multitask (precision)
            ppv_list, tpr_list = ef.compute_ppv_sensitivity(causal_fname,
                                                            selected,
                                                            args.num_features)
            with open(ppv_fname, 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in ppv_list]))
                f.close()

            with open(tpr_fname, 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in tpr_list]))
                f.close()
            #-----------------------------------------------------------
            #??? Does f.close() is necessary when usien with open statement ??? 

            #------------------------------------------------------------------
            # For each algorithm, for each task,
            # predict on the test set using a ridge-
            # regression trained with the selected features only.
            tr_indices = evalf.xp_indices[fold_idx]['trIndices']
            te_indices = evalf.xp_indices[fold_idx]['teIndices']

            for task_idx in range(args.num_tasks):
                # Single task
                fname = '%s/%s.sfan.fold_%d.task_%d.predicted' % \
                        (resu_dir, simu_id, fold_idx, task_idx)
                ef.run_ridge_selected(selected_st[task_idx], genotype_fname,
                                      phenotype_fname[task_idx],
                                      trIndices, teIndices, fname)

                # Multitask (no precision)
                fname = '%s/%s.msfan_np.fold_%d.task_%d.predicted' % \
                        (resu_dir, simu_id, fold_idx, task_idx)
                ef.run_ridge_selected(selected_np[task_idx], genotype_fname,
                                      phenotype_fname[task_idx],
                                      trIndices, teIndices, fname)

                # Multitask (precision)
                fname = '%s/%s.msfan.fold_%d.task_%d.predicted' % \
                        (resu_dir, simu_id, fold_idx, task_idx)
                ef.run_ridge_selected(selected[task_idx], genotype_fname,
                                      phenotype_fname[task_idx],
                                      trIndices, teIndices, fname)
            #------------------------------------------------------------------
        # END for fold_idx in range(args.num_folds)

                            
        #----------------------------------------------------------------------
        # TODO: For each algorithm, and for each task, compute RMSE
        # (define an external function, a bit as for ppv/tpr; use sklearn.metrics)
        # use the predictions saved in files (fold per fold)
        # the true values are given by phenotype_fname[task_idx] and te_indices
        # save to file '%s/%s.<algo>.rmse' % (args.resu_dir, simu_id)

        # Single task


        # Multitask (no precision)


        # Multitask (precision)                
        #----------------------------------------------------------------------




        #-----------------------------------------------------------------------
        # TODO: For each algorithm, and for each task, compute consistency index
        # between the features selected for each fold.
        # Use ef.consistency_index_k
        # (define an external function, a bit as above for RMSE)
        # use the selected features saved to files and the true causal features
        # save to file '%s/%s.<algo>.consistency' % (args.resu_dir, simu_id)

        # Single task


        # Multitask (no precision)


        # Multitask (precision)
        #-----------------------------------------------------------------------
    # END for repeat_idx in range(args.num_repeats)


    # TODO: Add line breaks in PPV and sensitivity files.


                            
    #--------------------------------------------------------------
    # TODO: For each algorithm compute average/mean for RMSE, PPVs,
    # sensitivities, consistency index.


    # TODO: Print out and save (with LaTeX table format) in
    # plain text file
    fname = '%s/%s.results' % (args.resu_dir, simu_id)
    #--------------------------------------------------------------

    

    #-------------
    # TODO: Plots
                            

    #-------------


if __name__ == "__main__":
    main()
        
