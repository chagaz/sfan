#!/bin/sh 
#$ -N repeats
python synthetic_data_experiments__parallel-repeat.py -k $1 -m $2 -n $3 -r $4 -f $5 -s $6 $7 $8 $9 -ri ${SGE_TASK_ID} --verbose
#python -c "from synthetic_data_experiment import run_repeat ; run_repeat()"
