import synthetic_data_experiments as sde
import argparse


if __name__ == "__main__":
    args = sde.get_integrous_arguments_values()
    analysis_files = sde.get_analysis_files_names(args.resu_dir, args.simu_id)
    sde.handle_measures_results(analysis_files, args)

