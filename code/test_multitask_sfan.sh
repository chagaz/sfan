clear

# Single task 
python multitask_sfan.py --num_tasks 1 \
                         --networks ../data/simu_multitask_01.network.dimacs \
                         --node_weights ../data/simu_multitask_01.scores_0.txt \
			 -l 0.001 -e 0.02 \
			 --output ../results/results_st_01.runtime.txt \
       > ../results/results_st_01.txt && cat ../results/results_st_01.txt && echo && echo

# Two tasks, but no correlation structure is given
python multitask_sfan.py --num_tasks 2 \
                         --networks ../data/simu_multitask_01.network.dimacs \
                         --node_weights ../data/simu_multitask_01.scores_0.txt \
                                        ../data/simu_multitask_01.scores_1.txt \
			 -l 0.001 -e 0.02 -m 0.01 \
			 --output ../results/results_mt.runtime_01.txt \
       > ../results/results_mt_01.txt && cat ../results/results_mt_01.txt && echo && echo

# Two tasks, a single network, given correlation structure
python multitask_sfan.py --num_tasks 2 \
                         --networks ../data/simu_multitask_01.network.dimacs \
                         --node_weights ../data/simu_multitask_01.scores_0.txt \
                                        ../data/simu_multitask_01.scores_1.txt \
                         --correlation_matrix ../data/simu_multitask_01.task_similarities.txt \
			 -l 0.001 -e 0.02 -m 0.01 \
			 --output ../results/results_mtc.runtime_01.txt \
       > ../results/results_mtc_01.txt && cat ../results/results_mtc_01.txt && echo && echo


# Two tasks, two networks, given correlation structure
python multitask_sfan.py --num_tasks 2 \
                         --networks ../data/simu_multitask_01.network.dimacs \
                                    ../data/simu_multitask_01.network.dimacs \
                         --node_weights ../data/simu_multitask_01.scores_0.txt \
                                        ../data/simu_multitask_01.scores_1.txt \
                         --correlation_matrix ../data/simu_multitask_01.task_similarities.txt \
			 -l 0.001 -e 0.02 -m 0.01 \
			 --output ../results/results_mtc.runtime_01.txt \
       > ../results/results_mtc_01.txt && cat ../results/results_mtc_01.txt && echo && echo

