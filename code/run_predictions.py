import synthetic_data_experiments as sde

import logging

if __name__ == "__main__":
    args = sde.get_integrous_arguments_values()
    
    for repeat_idx in xrange(args.num_repeats) : 
        
        resu_dir = "%s/repeat_%d" % (args.resu_dir, repeat_idx)
        data_dir = '%s/repeat_%d' % (args.data_dir, repeat_idx)
        
        trIndices_fname = data_dir+'/'+args.simu_id+'.fold%d.trIndices'
        teIndices_fname = data_dir+'/'+args.simu_id+'.fold%d.teIndices'
        ssIndices_fname = data_dir+'/'+args.simu_id+'.fold%d.ss%d.ssIndices'
        
        for fold_idx in xrange (args.num_folds) : 

            with open(trIndices_fname %(fold_idx), 'r') as trIndices_f : 
                line = trIndices_f.readline().split()
                trIndices = [int (i) for i in line ]
            with open(teIndices_fname %(fold_idx),'r') as teIndices_f : 
                line = teIndices_f.readline().split()
                teIndices =  [int (i) for i in line ]
                
            
            sde.run_predictions(fold_idx, args, resu_dir, data_dir, trIndices, teIndices )
    
