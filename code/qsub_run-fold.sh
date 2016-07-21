#!/bin/sh 
echo python synthetic_data_experiments__parallel-fold.py -k $1 -m $2 -n $3 -r $4 -f $5 -s $6 $7 $8 $9 $10 $11 $12 $SGE_TASK_ID --verbose  \
> $8'/SGE-outputs/snp'$2'r'$12'f'$SGE_TASK_ID.o 2>$8'/SGE-outputs/snp'$2'r'$12'f'$SGE_TASK_ID.e
python synthetic_data_experiments__parallel-fold.py -k $1 -m $2 -n $3 -r $4 -f $5 -s $6 $7 $8 $9 $10 $11 $12 $SGE_TASK_ID --verbose  \
> $8'/SGE-outputs/snp'$2'r'$12'f'$SGE_TASK_ID.o 2>$8'/SGE-outputs/snp'$2'r'$12'f'$SGE_TASK_ID.e
# "-k = num_tasks",
# "-m = num_features"
#"-n = num_samples"
# -r = num_repeats
#"-f = num_folds"
#"-s = num_subsamples"
#"data_dir"
#"resu_dir"
#"simu_id"
#"--hyperparam_fname_np"
#"--hyperparam_fname"
# repeat_idx"
# fold_idx"

