#!/bin/sh

#same code than run snp experimentation (experiment_numSNP.sh)

num_tasks=3
#num_features #range
num_samples=1000
num_repeats=10
num_folds=10
num_subsamples=20


#path='..' # Version BOB 
path=$SHAREDAT"/exp_changeNbSNP" #"/share/data40T/athenais/" # Version Cluster 


for num_features in 200 1000
do
    echo $num_features
    simu_id='numSNP_'$num_features
    data_dir=$path"/data/"$simu_id
    resu_dir=$path"/results_timing/"$simu_id
    
    python synthetic_data_experiments.py \
    -k $num_tasks -m $num_features -n $num_samples -r $num_repeats -f $num_folds -s $num_subsamples \
    $data_dir $resu_dir $simu_id --verbose > output$num_features

done

