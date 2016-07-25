import synthetic_data_experiments as sde
import evaluation_framework as ef
import logging
import plot
import numpy as np



def print_and_save_measure_table(measure_name, data, out_fname): 
    """ Print measures tables and save them (in LaTeX table format) in file

    Parameters
    ----------
    measure_name : string 
        measure name
    data : dict of list of list
        data[algo_name][task_idx][sample_idx] 
        = the value of the data 
        for the sample sample_idx, 
        for the task task_idx, 
        for the algo algo_name.
    out_fname : filename
        Path to the output file were LaTeX table has to be printed

    Side effects
    ------------
    Create a file named <out_fname> holding the table in LaTeX format.

    """

    header_print = (
            "-----------------------------------------------------------------------\n"
            "                      %s\n"
            "       +--------------------------------------------------------------+\n"
            "       |                              algo                            |\n"
            "  task |         SConES     |    MSConESnp       |        MSConES     |\n"
            "=======+====================+====================+====================+\n"
    )
    header_save = "task & sfan &np &msfan \\\\\\hline\n"
    algos_names = ('SConES', 'MSConESnp', 'MSConES')
    
    to_print = ''
    to_save = ''

    for task_idx in xrange (len (data[algos_names[0]] ) ) : 
        to_print += '{:^7d}|'.format(task_idx) 
        to_save += '{:^7d}&'.format(task_idx) 
        for algo in algos_names : 
            mean_task_algo = np.mean(data[algo][task_idx])
            std_task_algo = np.std(data[algo][task_idx])
            to_print += '{:9.3f} ±{:9.3f}|'.format(mean_task_algo, std_task_algo) 
            to_save += '{:9.3f} ±{:9.3f}&'.format(mean_task_algo, std_task_algo) 
        to_save = to_save[:-1] #don't want the last '&'
        to_save += "\\\\\n" # two step and not to_save[-1] = "\\\\\n" because python strings are immuable
        to_print+="\n"
        
    print  header_print  %measure_name + to_print

    with open(out_fname, 'w') as f : 
        f.write(header_save+to_save)




def extract_plotable_data_from_analysis_files(f_names, num_tasks, num_repeats, num_folds):
    """ Extract plotable data from analsis files 

    Arguments
    ---------
    f_names : filenames
        Path to files in the order st, np, msfan
        holding values whith space separated lists of values :
        -    one per task for each fold,
             or just one per task
        - one line per repeat
    num_tasks: int
        Number of tasks. 
    num_repeat : int 
        Number of repeat
    num_folds : int
        Number of folds.
        Needed for ppv and tpr files 
        because values per repeat are per fold then per task. 

    Return
    -------
    data : dict of list of list
        data[algo_name][task_idx][sample_idx] 
        = the value of the data 
        for the sample sample_idx, 
        for the task task_idx, 
        for the algo algo_name.
    """
    data = {}
    algos = ('SConES', 'MSConESnp', 'MSConES')

    # for acc, mcc, ppv and tpr file for which values per repeat = per line 
    # are per fold, then per task, 
    # we have to compute means per task taking account of value per repeat and per fold. 
    for algo_idx, algo in enumerate(algos) : 
        data_for_an_algo  = [[float() for i in xrange(num_repeats)] for j in xrange(num_tasks)]
        with open (f_names[algo_idx], 'r') as f : 
            # for repeat_idx, line in enumerate(f) :
            for repeat_idx in xrange(num_repeats) : 
                line = f.readline()
                line_content = [float (item) for item in line.split()]
                # line_content contains every float of the line 
                for task_idx in xrange(num_tasks) :
                    # we get values for the current task in content_for_a_task using a slice : 
                    # - if data are holded one line per repeat, then per task, for each fold, 
                    #   the slice take values stepped by num_tasks until the end of the line
                    # - if data are holded one line per repeat, one column per task as mean for every fold, 
                    #   then it works too, taking the task_idx th value
                    content_for_a_task = line_content[task_idx::num_tasks]
                    #print content_for_a_task
                    # content_for_a_task_contains float of the current task for each folds of the current repeat
                    # we take the mean of them to have the mean of the repeat
                    data_for_an_algo[task_idx][repeat_idx]  = np.mean(content_for_a_task)
            data[ algos[algo_idx] ] = data_for_an_algo
    return data



#####################################################################################################
#####################################################################################################
if __name__ == "__main__":
    """ Handle synthetic_data_experiments.py outputs

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



    1. Results are saved in following files : 

        <resu_dir>/<simu_id>.table.acc.tex
        <resu_dir>/<simu_id>.table.mcc.tex
        <resu_dir>/<simu_id>.table.ppv.tex
        <resu_dir>/<simu_id>.table.tpr.tex
        <resu_dir>/<simu_id>.table.rmse.tex
        <resu_dir>/<simu_id>.table.ci.tex
            Average/standard deviation values for: consistency index, RMSE, PPV and TPR, as LaTeX table.



    For each algo in (`sfan`, `msfan_np`, `msfan`):

        <resu_dir>/<simu_id>.<algo>.rmse : 
          List of final RMSEs 
          one line per repeat
          for each repeat, one value per task
        <resu_dir>/<simu_id>.<algo>.consistency : 
          List of final Consistency Indices 
          one line per repeat
          for each repeat, one value per task

        For each classification measure 
        (accuracy (acc), Mathieu coefficient (mcc), Prositive Predictive Value (ppv) and True Positive Value (tpr) )
        a file named `<resu_dir>/<simu_id>.<algo>.<measure>` :
        <args.resu_dir>/<simu_id>.<algo>.acc
        <args.resu_dir>/<simu_id>.<algo>.mcc
        <args.resu_dir>/<simu_id>.<algo>.ppv
        <args.resu_dir>/<simu_id>.<algo>.tpr 
            Space-separated lists of PPVs (one value per task and per fold),
            each line corresponds to one repeat. 


        #### For each repeat, fold and task : 
        <resu_dir>/<repeat_idx>/<simu_id>.<algo>.fold_<fold_idx>.task_<task_idx>.predicted : 
            phenotype prediction of the test set using a ridge-regression trained with the selected features only.
            One value per line, one line per sample
    
    2. Charts :
        For each measure (ci, mcc, ppv, tpr, rmse, acc), 
        a chart named <simu_id>.<measure>.png : one boxplot per algorithm, 
        grouped per task, with error bars. 
        <args.resu_dir>/<simu_id>.ci.png
        <args.resu_dir>/<simu_id>.rmse.png
        <args.resu_dir>/<simu_id>.acc.png
        <args.resu_dir>/<simu_id>.mcc.png
        <args.resu_dir>/<simu_id>.ppv.png
        <args.resu_dir>/<simu_id>.tpr.png


    For file format specifications, see README.md


    Example
    -------
    $ python handle-output.py -k 3 -m 200 -n 100 -r 10 -f 10 -s 10 \
             ../data/simu_synth_01 ../results/simu_synth_01 simu_01 --verbose


    """
    args = sde.get_integrous_arguments_values()

    for repeat_idx in xrange(args.num_repeats) : 
        print '=========== repeat : ', repeat_idx
        resu_dir = "%s/repeat_%d" % (args.resu_dir, repeat_idx)
        data_dir = '%s/repeat_%d' % (args.data_dir, repeat_idx)
        

        analysis_files = sde.get_analysis_files_names(args.resu_dir, args.simu_id)
        #-----------------
        genotype_fname = '%s/%s.genotypes.txt' % (data_dir, args.simu_id)
        #------------------
        phenotypes_fnames = ['%s/%s.phenotype_%d.txt' % \
                        (data_dir, args.simu_id, task_idx) \
                        for task_idx in range(args.num_tasks)]

        #------------------
        # get causal features from files
        causal_fname = '%s/%s.causal_features.txt' % (data_dir, args.simu_id)
        causal_features = []
        with open(causal_fname, 'r') as f:
            for line_idx, line in enumerate(f):
                causal_features.append( map(int, line.split()) )
        #-------------------
        

        #-------------------
        trIndices_fname = resu_dir+'/'+args.simu_id+'.fold%d.trIndices'
        teIndices_fname = resu_dir+'/'+args.simu_id+'.fold%d.teIndices'
        ssIndices_fname = resu_dir+'/'+args.simu_id+'.fold%d.ss%d.ssIndices'
        
        xp_indices = [{'trIndices': list(), 'teIndices':list(), 'ssIndices':list()} for fold in xrange(args.num_folds)]

        for fold_idx in xrange (args.num_folds) : 
            print 'folds : ', fold_idx
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
            # print in analysis files 
            # Files structure : 
            # 1 line per repeat
            # on each line : valTask1, valTask2, ... valTaskn for each fold
            
            with open(analysis_files['acc_st'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in acc_list_st]))
            with open(analysis_files['acc_msfan_np'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in acc_list_np]))
            with open(analysis_files['acc_msfan'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in acc_list_msfan]))

            with open(analysis_files['mcc_st'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in mcc_list_st]))
            with open(analysis_files['mcc_msfan_np'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in mcc_list_np]))
            with open(analysis_files['mcc_msfan'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in mcc_list_msfan]))

            with open(analysis_files['ppv_st'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in ppv_list_st]))
            with open(analysis_files['ppv_msfan_np'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in ppv_list_np]))
            with open(analysis_files['ppv_msfan'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in ppv_list_msfan]))

            with open(analysis_files['tpr_st'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in tpr_list_st]))
            with open(analysis_files['tpr_msfan_np'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in tpr_list_np]))
            with open(analysis_files['tpr_msfan'], 'a') as f:
                f.write('%s ' % ' '.join(['%.2f ' % x for x in tpr_list_msfan]))
                
            #-----------------------------------------------------------------------   
            # Run predictions : 
             
            # For each algorithm, for each task,
            # predict on the test set using a ridge-
            # regression trained with the selected features only.

            for task_idx in xrange(args.num_tasks):
                print 'task', task_idx
                logging.info ('task n. %d' %task_idx)
                # Single task
                logging.info("prediction st")
                fname = '%s/%s.sfan.fold_%d.task_%d.predicted' % \
                        (resu_dir, args.simu_id, fold_idx, task_idx)
                ef.run_ridge_selected(selected_st[task_idx], genotype_fname,
                                      phenotypes_fnames[task_idx],
                                      xp_indices[fold_idx]['trIndices'], xp_indices[fold_idx]['teIndices'], fname)

                # Multitask (no precision)
                logging.info("pred np")
                fname = '%s/%s.msfan_np.fold_%d.task_%d.predicted' % \
                        (resu_dir, args.simu_id, fold_idx, task_idx)
                ef.run_ridge_selected(selected_np[task_idx], genotype_fname,
                                      phenotypes_fnames[task_idx],
                                      xp_indices[fold_idx]['trIndices'], xp_indices[fold_idx]['teIndices'], fname)

                # Multitask (precision)
                logging.info("pred msfan")
                fname = '%s/%s.msfan.fold_%d.task_%d.predicted' % \
                        (resu_dir, args.simu_id, fold_idx, task_idx)
                ef.run_ridge_selected(selected[task_idx], genotype_fname,
                                      phenotypes_fnames[task_idx],
                                      xp_indices[fold_idx]['trIndices'], xp_indices[fold_idx]['teIndices'], fname)

        # END for fold_idx in range(args.num_folds)
        
        #----------------------------------------------------------------------
        logging.info( "======== Compute & save RMSE")
        # For each algorithm, and for each task, compute RMSE
        # using :
        #   - an external function : compute_ridge_selected_RMSE() returns list of rmse, task per task
        #   - the predictions saved in files (fold per fold)
        #   - the true values given by phenotypes_fnames[task_idx] and te_indices
        # save to file '%s/%s.<algo>.rmse' % (args.resu_dir, args.simu_id)
        # => rmse_st_fname ; rmse_np_fname ; rmse_fname
        # Files structure : 
        # each line = a repeat
        # on each line there are several RMSE values, one per task

        # Single task
        predicted_phenotypes_fname = resu_dir+'/'+args.simu_id+'.sfan.fold_%d.task_%d.predicted' 
        print predicted_phenotypes_fname
        print xp_indices
        rmse_list = ef.compute_ridge_selected_RMSE( phenotypes_fnames, predicted_phenotypes_fname, 
                                        xp_indices, args.num_tasks)
        print rmse_list
        with open(analysis_files['rmse_st'], 'a') as f:
            f.write('%s \n' % ' '.join(['%.2f ' % x for x in rmse_list]))
        # Multitask (no precision)
        predicted_phenotypes_fname = resu_dir+'/'+args.simu_id+'.msfan_np.fold_%d.task_%d.predicted' 
        rmse_list = ef.compute_ridge_selected_RMSE( phenotypes_fnames, predicted_phenotypes_fname, 
                                        xp_indices, args.num_tasks)
        print rmse_list
        with open(analysis_files['rmse_msfan_np'], 'a') as f:
            f.write('%s \n' % ' '.join(['%.2f ' % x for x in rmse_list]))
        # Multitask (precision)
        predicted_phenotypes_fname = resu_dir+'/'+args.simu_id+'.msfan.fold_%d.task_%d.predicted' 
        rmse_list = ef.compute_ridge_selected_RMSE( phenotypes_fnames, predicted_phenotypes_fname, 
                                        xp_indices, args.num_tasks)       
        print rmse_list      
        with open(analysis_files['rmse_msfan'], 'a') as f:
            f.write('%s \n' % ' '.join(['%.2f ' % x for x in rmse_list]))
        #----------------------------------------------------------------------

        #-----------------------------------------------------------------------
        logging.info( "======== Compute & save CI ")
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
        # At the end of a repeat, add a carriage return to analysis files holding measures assessing classification : 
        # acc, mcc, ppv, tpr files : 
        with open(analysis_files['acc_st'], 'a') as f:
            f.write('\n')
        with open(analysis_files['acc_msfan_np'], 'a') as f:
            f.write('\n')
        with open(analysis_files['acc_msfan'], 'a') as f:
            f.write('\n')

        with open(analysis_files['mcc_st'], 'a') as f:
            f.write('\n')
        with open(analysis_files['mcc_msfan_np'], 'a') as f:
            f.write('\n')
        with open(analysis_files['mcc_msfan'], 'a') as f:
            f.write('\n')

        with open(analysis_files['ppv_st'], 'a') as f:
            f.write('\n')
        with open(analysis_files['ppv_msfan_np'], 'a') as f:
            f.write('\n')
        with open(analysis_files['ppv_msfan'], 'a') as f:
            f.write('\n')

        with open(analysis_files['tpr_st'], 'a') as f:
            f.write('\n')
        with open(analysis_files['tpr_msfan_np'], 'a') as f:
            f.write('\n')
        with open(analysis_files['tpr_msfan'], 'a') as f:
            f.write('\n')

        
    # END for repeat_idx in xrange(args.num_repeats)
    
    
    #----------------------------------------------------------------------------
    # For each measures :
    #  - extract data of measures from analysis_files, 
    #  - print tables, 
    #  - save tables in Latex format, 
    #  - and Plots
    
    template_name =  args.resu_dir+'/'+args.simu_id+'.%s.png'
    
    # data [algo_id][task_id] = list of values

    data = extract_plotable_data_from_analysis_files(
        [analysis_files['acc_st'], analysis_files['acc_msfan_np'], analysis_files['acc_msfan'] ], 
        args.num_tasks, args.num_repeats, args.num_folds
    )
    print_and_save_measure_table("accuracy", data, "%s/%s.table.accuracy.tex"%(args.resu_dir, args.simu_id ))
    plot.horizontal_boxplots(data, template_name%'accuracy')

    data = extract_plotable_data_from_analysis_files(
        [analysis_files['mcc_st'], analysis_files['mcc_msfan_np'], analysis_files['mcc_msfan'] ], 
        args.num_tasks, args.num_repeats, args.num_folds
    )
    print_and_save_measure_table("mcc", data, "%s/%s.table.mcc.tex"%(args.resu_dir, args.simu_id ))
    plot.horizontal_boxplots(data, template_name%'mcc')

    data = extract_plotable_data_from_analysis_files(
        [analysis_files['ppv_st'], analysis_files['ppv_msfan_np'], analysis_files['ppv_msfan'] ], 
        args.num_tasks, args.num_repeats, args.num_folds
    )
    print_and_save_measure_table("ppv", data, "%s/%s.table.ppv.tex"%(args.resu_dir, args.simu_id ))
    plot.horizontal_boxplots(data, template_name%'ppv')

    data = extract_plotable_data_from_analysis_files(
        [analysis_files['tpr_st'], analysis_files['tpr_msfan_np'], analysis_files['tpr_msfan'] ], 
        args.num_tasks, args.num_repeats, args.num_folds
    )
    print_and_save_measure_table("tpr", data, "%s/%s.table.tpr.tex"%(args.resu_dir, args.simu_id ))
    plot.horizontal_boxplots(data, template_name%'tpr')

    data = extract_plotable_data_from_analysis_files(
        [analysis_files['rmse_st'], analysis_files['rmse_msfan_np'], analysis_files['rmse_msfan'] ], 
        args.num_tasks, args.num_repeats, args.num_folds
    )
    print_and_save_measure_table("rmse", data, "%s/%s.table.rmse.tex" %(args.resu_dir, args.simu_id ))
    plot.horizontal_boxplots(data, template_name%'rmse')

    data = extract_plotable_data_from_analysis_files(
        [analysis_files['ci_st'], analysis_files['ci_msfan_np'], analysis_files['ci_msfan'] ], 
        args.num_tasks, args.num_repeats, args.num_folds
    )
    print_and_save_measure_table("ci", data,  "%s/%s.table.ci.tex" %(args.resu_dir, args.simu_id ))
    plot.horizontal_boxplots(data, template_name%'ci')


