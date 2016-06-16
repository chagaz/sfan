import synthetic_data_experiments as sde
import logging

if __name__ == "__main__":
    args = sde.get_integrous_arguments_values()

    for repeat_idx in xrange(args.num_repeats) : 
        
        resu_dir = "%s/repeat_%d" % (args.resu_dir, repeat_idx)
        data_dir = '%s/repeat_%d' % (args.data_dir, repeat_idx)
        
        #------------------
        # get causal features from files
        causal_fname = '%s/%s.causal_features.txt' % (data_dir, args.simu_id)
        causal_features = []
        with open(causal_fname, 'r') as f:
            for line_idx, line in enumerate(f):
                causal_features.append( map(int, line.split()) )
        #-------------------
        

        #-------------------
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
            
            
            #-----------------------------------------------------------------------------
            # Get selected features from files : 
            
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

            #--------------------------------------------------------------------------------
            # Measure computation : 
            
            # For each algorithm, and for each task, compute measures
            #   - accuracy_score
            #   - matthews_corrcoef
            #   - precision_score = pvv
            #   - recall_score = tpr
            # For the current repeat and the current fold, 
            # ppv_list ant tpr_list and list of ppv and tpr respectively
            # for each task
            

            # Single task
            acc_list_st, mcc_list_st, ppv_list_st, tpr_list_st = ef.evaluate_classification(
                                                            causal_features,
                                                            selected_st,
                                                            args.num_features)
            # Multitask (no precision)
            acc_list_np, mcc_list_np, ppv_list_np, tpr_list_np = ef.evaluate_classification(
                                                            causal_features,
                                                            selected_np,
                                                            args.num_features)
            # Multitask (precision)
            acc_list_msfan, mcc_list_msfan, ppv_list_msfan, tpr_list_msfan = ef.evaluate_classification(
                                                            causal_features,
                                                            selected,
                                                            args.num_features)
            #--------------------------------------------------------------------------------
            # Measure saving in <measure>_fname
                                                #      v-- = fold_id
            # print in files : %s/repeat%d/%s.sfan.fold*.ppv" % (resu_dir, repeat_idx, simu_id)
            # with '%.2f' precision
            # Files structure : 
            # 1 line per repeat
            # on each line : valTask1, valTask2, ... valTaskn for each fold
            acc_template_f_name = str(resu_dir)+"/"+str(args.simu_id)+".%s.fold_"+str(fold_idx)+".acc"
            mcc_template_f_name = str(resu_dir)+"/"+str(args.simu_id)+".%s.fold_"+str(fold_idx)+".mcc"
            ppv_template_f_name = str(resu_dir)+"/"+str(args.simu_id)+".%s.fold_"+str(fold_idx)+".ppv"
            tpr_template_f_name = str(resu_dir)+"/"+str(args.simu_id)+".%s.fold_"+str(fold_idx)+".tpr"

            # Single task

            with open (acc_template_f_name %"sfan", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in acc_list_st]))
            with open (mcc_template_f_name %"sfan", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in mcc_list_st]))
            with open (ppv_template_f_name %"sfan", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in ppv_list_st]))
            with open (tpr_template_f_name %"sfan", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in ppv_list_st]))

            # Multitask (no precision)

            with open (acc_template_f_name %"msfan_np", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in acc_list_st]))
            with open (mcc_template_f_name %"msfan_np", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in mcc_list_st]))
            with open (ppv_template_f_name %"msfan_np", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in ppv_list_np]))
            with open (tpr_template_f_name %"msfan_np", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in tpr_list_np]))

            # Multitask (precision)

            with open (acc_template_f_name %"msfan", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in acc_list_st]))
            with open (mcc_template_f_name %"msfan", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in mcc_list_st]))
            with open (ppv_template_f_name %"msfan", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in ppv_list_msfan]))
            with open (tpr_template_f_name %"msfan", 'w') as f:
                f.write('%s \n' % ' '.join(['%.2f ' % x for x in tpr_list_msfan]))
                
            #-----------------------------------------------------------------------   
            # Run predictions : 
             
            # For each algorithm, for each task,
            # predict on the test set using a ridge-
            # regression trained with the selected features only.

            for task_idx in xrange(args.num_tasks):
                logging.info ('task n. %d' %task_idx)
                # Single task
                logging.info("prediction st")
                fname = '%s/%s.sfan.fold_%d.task_%d.predicted' % \
                        (resu_dir, args.simu_id, fold_idx, task_idx)
                ef.run_ridge_selected(selected_st[task_idx], genotype_fname,
                                      phenotype_fnames[task_idx],
                                      trIndices, teIndices, fname)

                # Multitask (no precision)
                logging.info("pred np")
                fname = '%s/%s.msfan_np.fold_%d.task_%d.predicted' % \
                        (resu_dir, args.simu_id, fold_idx, task_idx)
                ef.run_ridge_selected(selected_np[task_idx], genotype_fname,
                                      phenotype_fnames[task_idx],
                                      trIndices, teIndices, fname)

                # Multitask (precision)
                logging.info("pred msfan")
                fname = '%s/%s.msfan.fold_%d.task_%d.predicted' % \
                        (resu_dir, args.simu_id, fold_idx, task_idx)
                ef.run_ridge_selected(selected[task_idx], genotype_fname,
                                      phenotype_fnames[task_idx],
                                      trIndices, teIndices, fname)

        # END for fold_idx in range(args.num_folds)
        
        sde.print_analysis_files( args, resu_dir, data_dir, xp_indices)
        # concatenation ppv, 
        # concatenation tpr
        # compute & write RMSE 
        # compute & write CI
        
    # END for repeat_idx in xrange(args.num_repeats)
    
    
    #----------------------------------------------------------------------------
    
    analysis_files = sde.get_analysis_files_names(args.resu_dir, args.simu_id)
    
    # For each measure compute average/mean +- standard deviation per task 
    means, std = extract_res_means_and_std(analysis_files, args)
    
    
    # Write means and std : 
    fname = args.resu_dir+'/'+args.simu_id+'.results_%s'
    print_save_res_measures(means, std, fname)

    # Plots : 

    for measure in means : 
        f_name = "%s/%s.%s_plot.values" %(args.resu_dir, args.simu_id, measure)
        print_plot_files(f_name, means[measure], std[measure])
        plot.bar_plot(measure, f_name)
    
    
