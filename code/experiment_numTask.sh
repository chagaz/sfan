#!/bin/sh



###############################################################################
# var to modify : 

# -- param of the experimentation : 
array_tasks=( 3 4 5 6 7 8 9 10 11 12 ) # range of num_tasks values
num_features=200
num_samples=1000
num_repeats=10
num_folds=10
num_subsamples=20

# -- about the paths : 

#path='..' # Version BOB 
path="/share/data40T/athenais/" # Version Cluster 

dat_path=$path'dat/experiment_numTask/'
res_path=$path'res/experiment_numTask_'`date +%Y-%m-%d-%H:%M`'/'


###############################################################################



for num_tasks in "${array_tasks[@]}"
do
    echo $num_tasks
    simu_id='numTask_'$num_tasks
    
    data_dir=$dat_path'/'$simu_id'/'
    resu_dir=$res_path'/'$simu_id'/'
    
    python generate_data.py -k $num_tasks  -m $num_features -n $num_samples $data_dir $simu_id
    
    python synthetic_data_experiments.py \
    -k $num_tasks -m $num_features -n $num_samples -r $num_repeats -f $num_folds -s $num_subsamples \
    $data_dir $resu_dir $simu_id

done

