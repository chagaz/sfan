import synthetic_data_experiments as sde
import argparse
import logging

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
    
    parser.add_argument("hyperparam_fname", help="Number of CV folds")
                                  # arg that differ with sde. 
    parser.add_argument("repeat_idx", help="Index of the current repeat",
                        type=int) # arg that differ with sde. 
    parser.add_argument("fold_idx", help="Index of the current fold",
                        type=int) # arg that differ with sde. 

    parser.add_argument("-v", "--verbose", help="Turn on detailed info log",
                        action='store_true')

    args = parser.parse_args()
    
    args.fold_idx = args.fold_idx -1
    
    if args.verbose:
        logging.basicConfig(format="[%(levelname)s] %(message)s",
                            level=logging.DEBUG)
        logging.info("Verbose output.")


    resu_dir = "%s/repeat_%d" % (args.resu_dir, args.repeat_idx)
    
    data_dir = '%s/repeat_%d' % (args.data_dir, args.repeat_idx)
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
    
    
    #with open(args.hyperparam_fname) as f:
        #lbd_eta_mu_values = f.readlines()
    #lbd_eta_values = [" ".join(plist.split()[:-2]) \
    #                  for plist in lbd_eta_mu_values]
    
    lbd_eta_mu_values = []
    lbd_eta_values = []
    with open(args.hyperparam_fname) as f:
        for line in f : 
            lbd_eta_mu_values.append(line)
            lbd_eta_values.append(" ".join(line.split()[:-2]) )
    
    
    # indices for this fold : 
    # TODO : factorisation of fname template...
    trIndices_fname = data_dir+'/'+args.simu_id+'.fold%d.trIndices'
    teIndices_fname = data_dir+'/'+args.simu_id+'.fold%d.teIndices'
    ssIndices_fname = data_dir+'/'+args.simu_id+'.fold%d.ss%d.ssIndices'

    indices = {'trIndices': list(), 'teIndices':list(), 'ssIndices':list()}
    with open(trIndices_fname %(args.fold_idx), 'r') as trIndices_f : 
        line = trIndices_f.readline().split()
        indices["trIndices"] = [int (i) for i in line ]
    with open(teIndices_fname %(args.fold_idx),'r') as teIndices_f : 
        line = teIndices_f.readline().split()
        indices["teIndices"] =  [int (i) for i in line ]
    for ss_idx in xrange(args.num_subsamples) :
        with open(ssIndices_fname  %(args.fold_idx,ss_idx), 'r') as ssIndices_f:
            line = ssIndices_f.readline().split()
            indices["ssIndices"].append(  [int (i) for i in line ] )        
                      
    tmp_weights_fnames = sde.get_tmp_weights_fnames(resu_dir, args.simu_id, args.fold_idx)
             
    sde.run_fold(
            args.fold_idx,
            args, 
            lbd_eta_values, lbd_eta_mu_values, 
            indices, 
            genotype_fname, network_fname , tmp_weights_fnames, precision_fname , causal_fname, phenotype_fnames, scores_fnames,
            resu_dir)
    
