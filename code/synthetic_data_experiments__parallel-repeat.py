import synthetic_data_experiments as sde
import argparse


if __name__ == "__main__":
    
    # TODO : use sde.get_integrous_arg_values ???
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
                        
    parser.add_argument("-ri", "--repeat_idx", help="Index of the current repeat",
                        type=int) #the only arg that differ with sde. 



    args = parser.parse_args()
    analysis_files = sde.get_analysis_files_names(args.resu_dir, args.simu_id)
    sde.run_repeat(args.repeat_idx, args, analysis_files)
    
