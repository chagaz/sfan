#!/bin/sh
# inspired of http://ccn.ucla.edu/wiki/index.php/How_to_have_a_script_wait_for_jobs_before_proceeding


sleep_time=900 # seconds; don't make this too short! don't want to tax system with excessive qstat calls
# 900 seconds = 15 mins


num_tasks=3
#num_features #range
num_samples=1000
num_repeats=10
num_folds=10
num_subsamples=20


#path='..' # Version BOB 
path="/share/data40T/athenais/exp_numSNP_30juin" # Version Cluster 


for num_features in 200 1000 10000 20000 50000
do
    echo $num_features
    simu_id='numSNP_'$num_features
    data_dir=$path"/data/"$simu_id
    resu_dir=$path"/results/"$simu_id

    # hold before running another python script...
    counter=`qstat | grep " r "|wc -l` # how many jobs are currently running
    until [ $counter -ge 100 ] # while number of running job not <= 100
    do
	    sleep $sleep_time
	    counter=`qstat | grep " r "|wc -l`
    done
	
	python synthetic_data_experiments.py \
    -k $num_tasks -m $num_features -n $num_samples -r $num_repeats -f $num_folds -s $num_subsamples \
    $data_dir $resu_dir $simu_id --verbose


done

