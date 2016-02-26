"""synthetic_data_experiments.py -- Run validation experiments on synthetic data

In this version all experiments are run sequentially.
"""

import argparse
import logging
import os
import subprocess
import sys

import generate_data
import multitask_sfan


def main():
    """
    Sequentially run validation experiments on synthetic data.

    Arguments
    ---------


    Outputs
    -------


    Example
    -------
    $ python synthetic_data_experiments -v -k 2 -m 1000 -n 150 -r 10 -f 10 -b 10 \
             ../data/simu_synth_01 ../results/simu_synth_01 simu_01
    """
    # Get arguments values
    parser = argparse.ArgumentParser(description="Validation experiments on synthetic data",
                                     add_help=True)
    parser.add_argument("-v", "--verbose", help="Turn on more detailed info log",
                        action='store_true')
    parser.add_argument("-k", "--num_tasks", help="Number of tasks", type=int)
    parser.add_argument("-m", "--num_features", help="Number of features", type=int)
    parser.add_argument("-n", "--num_samples", help="Number of samples", type=int)
    parser.add_argument("-r", "--num_repeats", help="Number of repeats", type=int)
    parser.add_argument("-f", "--num_folds", help="Number of CV folds", type=int)
    parser.add_argument("-b", "--num_bootstraps", help="Number of bootstrap iterations",
                        type=int)
    parser.add_argument("data_dir", help="Simulated data directory")
    parser.add_argument("res_dir", help="Results directory")
    parser.add_argument("simu_id", help="Simulation name")
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
        logging.error("The number of features must be larger than NUM_CAUSAL_TOTAL (%d).\n" % \
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
        assert(args.num_bootstraps > 0)
    except AssertionError:
        logging.error("The number of bootstraps must be strictly positive\n")
        logging.error("Use --help for help.\n")
        sys.exit(-1)

    # Verbose
    if args.verbose:
        logging.basicConfig(format="[%(levelname)s] %(message)s", level=logging.DEBUG)
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
    if not os.path.isdir(args.res_dir):
        logging.info("Creating %s\n" % args.res_dir)
        try: 
            os.makedirs(args.res_dir)
        except OSError:
            if not os.path.isdir(args.res_dir):
                raise


    for repeat_idx in range(args.num_repeats):
        # Instantiate data generator
        data_gen = generate_data.SyntheticDataGenerator(args.num_tasks, args.num_features,
                                                        args.num_samples,
                                                        args.data_dir, args.simu_id)
        # Generate modular data
        data_gen.generate_modular()

        # Name of data files
        # Either "hard-code" or edit generate_modular to return the names of these files
        network_fname = '%s/%s.network.dimacs' % (args.data_dir, args.simu_id)
        correlation_fname = TODO
        phenotype_fnames = TODO
        

        # Instantiate evaluation framework
        ef = evaluation_framework.Framework(args.num_samples, args.num_folds, args.num_bootstrapsxs)

        # Compute cross-validation folds and bootstrap indices
        ef.compute_indices()

        # Save cross-validation folds and bootstrap indices to file
        ef.save_indices(args.data_dir, args.simu_id)

        # TODO: define the ranges of values for lbd, eta, mu
        # Note: acceptable values of lbd, eta and mu are linked,
        # we won't define all three ranges as constants.

        # Note: MU_VALUES must contain 0, as this corresponds to
        # solving all tasks separately (single-task)
        lbd_values = []
        eta_values = []
        mu_values = [0]
        
        
        
        for fold_idx in range(args.num_folds):
            # Intitalize dictionary to store selected features
            # selected_features is a list of lists of selected features (one per bootstrap iteration)
            # selected_features_dict is a nested dictionary, indexed by
            #   - value of mu
            #   - value of lbd
            #   - value of eta
            #   - value of task_idx
            # i.e. you get a specific selected_features from it by querying
            # selected_features_dict[str(mu)][str(lbd)][str(eta)][task_idx]
            selected_features_dict = {}          
            selected_features_nocorr_dict = {}
            for lbd in range(lbd_values):
                selected_features_dict[str(lbd)] = {}          
                selected_features_nocorr_dict[str(lbd)] = {}
                for eta in range(eta_values):
                    selected_features_dict[str(lbd)][str(eta)] = {}          
                    selected_features_nocorr_dict[str(lbd)][str(eta)] = {}
                    for mu in range(mu_values):
                        selected_features_dict[str(lbd)][str(eta)][str(mu)] = {}          
                        selected_features_nocorr_dict[str(lbd)][str(eta)][str(mu)] = {}
                        for task_idx in range(args.num_tasks):
                            selected_features_dict[str(lbd)][str(eta)][str(mu)][task_idx] = []          
                            selected_features_nocorr_dict[str(lbd)][str(eta)][str(mu)][task_idx] = []
                            
                
            for bs_idx in range(args.num_bootstraps):
                # Get samples
                sample_indices = ef.xp_indices[fold_idx]['bsIndices'][bs_idx]

                # TODO:
                # Generate sample-specific network scores
                # Save to files
                weights_fname = []


                
                for lbd in range(lbd_values):
                    for eta in range(eta_values):
                        for mu in range(mu_values):

                            sfnd_ = selected_features_nocorr_dict[str(mu)][str(lbd)][str(eta)]
                            sfd_ = selected_features_dict[str(mu)][str(lbd)][str(eta)]
                            
                            # Run multi-task sfans *without* correlation matrix
                            # Ideally, I'd do the following:
                            # sfan_solver = Sfan(args.num_tasks, network_fname, weights_fname,
                            #                    lbd, eta, mu, correlation_fname, output_f='/tmp/test')
                            # tt = sfan_solver.create_dimacs()
                            # sfan_solver.run_maxflow()

                            # But because cython output to screen is NOT caught by sys.stdout,
                            # we need to run this externally -- it's a bit ugly.
                            p = subprocess.Popen(['python', 'multitask_sfan.py',
                                                  '--num_tasks', str(args.num_tasks),
                                                  '--networks', network_fname,
                                                  '--node_weights', weights_fnames,
                                                  '-l', str(lbd), '-e', str(eta), '-m', str(mu),
                                                  '--output', '/tmp/test'], # except if we're saving runtimes!
                                                  stdout=subprocess.PIPE)
                            selected_ = [[int(x) for x in line.split()] \
                                         for line in p.communicate()[0].split("\n")[2:]]

                            # Store selected_features_nocorr in the dictionary
                            for task_idx, selected_list in enumerate(selected_):
                                sfnd_[task_idx].append(selected_list)

                            # run multi-task sfans *with* correlation matrix
                            p = subprocess.Popen(['python', 'multitask_sfan.py',
                                                  '--num_tasks', str(args.num_tasks),
                                                  '--networks', network_fname,
                                                  '--node_weights', weights_fnames,
                                                  '--correlation_matrix', correlation_fname,
                                                  '-l', str(lbd), '-e', str(eta), '-m', str(mu),
                                                  '--output', '/tmp/test'], # except if we're saving runtimes!
                                                  stdout=subprocess.PIPE)
                            selected_ = [[int(x) for x in line.split()] \
                                         for line in p.communicate()[0].split("\n")[2:]]

                            # Store selected_features in the dictionary
                            for task_idx, selected_list in enumerate(selected_):
                                sfd_[task_idx].append(selected_list)
            # END for bs_idx in range(args.num_bootstraps)
                            
            # Compute consistencies for single-task, multi-task with no correlation matrix,
            # multi-task with correlation matrices. Get optimal parameter values for each algo.
            for lbd in range(lbd_values):
                for eta in range(eta_values):

                    # Single-task: mu = 0
                    for task_idx in range(args.num_tasks):
                        selected_list = selected_features_dict[str(0.)][str(lbd)][str(eta)][task_idx]
                        cidx = evaluation_framework.consistency_index_k(selected_list)
                        # TODO ...

                    # Average consistency indices across tasks
                    # TODO ...
                    
                    for mu in range(mu_values):
                        selected_list_d = selected_features_nocorr_dict[str(mu)]
                        for task_idx in range(args.num_tasks):
                            selected_list = selected_list_d[str(lbd)][str(eta)][task_idx]
                            cidx = evaluation_framework.consistency_index_k(selected_list)
                            # TODO ...

                    for mu in range(mu_values):
                        selected_list_d = selected_features_dict[str(mu)]
                            selected_list = selected_list_d[str(lbd)][str(eta)][task_idx]
                            cidx = evaluation_framework.consistency_index_k(selected_list)
                            # TODO ...

                    # Average consistency indices across tasks


            # Run algorithms again to select features, using the whole training set


            # Train a ridge-regression on the selected features / training set
            # Use sklearn.linear_models


            # Predict on test set



            # Compute PPV and sensitivity
            # Use sklearn.metrics

                            
        # END for fold_idx in range(args.num_folds)

    # Compute RMSE 


    # Compute CI



    # END for repeat_idx in range(args.num_repeats)


    # Save RMSE, PPVs, sensitivities, in plain text files

                            
    # Compute average/mean for these values                         
    # Print out and save in .tex tables


    # Display (use matplotlib)
                            


        