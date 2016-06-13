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
        
        xp_indices = [{'trIndices': list(), 'teIndices':list(), 'ssIndices':list()} for fold in xrange(args.num_folds)]

        for fold_idx in xrange (args.num_folds) : 
            #----------------------------------------------------------------------------
            # get xp_indices from files : 
            
            with open(trIndices_fname %(fold_idx), 'r') as trIndices_f : 
                line = trIndices_f.readline().split()
                xp_indices[fold_idx]["trIndices"] = [int (i) for i in line ]
            with open(teIndices_fname %(fold_idx),'r') as teIndices_f : 
                line = teIndices_f.readline().split()
                xp_indices[fold_idx]["teIndices"] =  [int (i) for i in line ]
            
            for ss_idx in xrange (args.num_subsamples) : 
                with open(ssIndices_fname  %(fold_idx,ss_idx), 'r') as ssIndices_f:
                    line = ssIndices_f.readline().split()
                    xp_indices[fold_idx]["ssIndices"].append( [int (i) for i in line ] ) 
            #----------------------------------------------------------------------------

