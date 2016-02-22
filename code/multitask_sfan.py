""" multitask_sfan.py: Generate the super-network for the given multi-task network-guided
feature selection problem and run maxflow on it.
"""

import argparse
import sys
import numpy as np
from numpy.linalg import inv
import time

import gt_maxflow



def main() : 
    """ Solve a multi-task network-guided feature selection problem by
    generating the corresponding super-network and runing maxflow on it.

    Arguments
    ---------
    args.num_tasks: int
        Number of tasks.
    args.networks: filename(s)
        List of paths of network files, space-separated.
        If the same network is to be used for all files, only specify one file.
    args.nodes_weight: filename(s)
        List of paths of network node weights (i.e. feature relevance scores), space-separated.
        One list per task.
    args.correlation_matrix: {filename, None}, optional
        Path of correlation matrix.
    args.lambdA: float
        Regularization parameter for connectivity.
    args.eta: float
        Regularization paramter for sparsity.
    args.mu: {float, None}, optional
        Regularization paramter for task relatednes.
    args.output: {filename, None}, optional
        File where to store computation run times.

    Output
    ------
    Printed to screen, the list of nodes selected in each network. 
    

    Examples
    --------
    Single task (equivalent to SConES [1]):
    $ python multitask_sfan.py --num_task 1 \
             --networks ../data/simu_multitask_01.network.dimacs \
             --nodes_weight ../data/simu_multitask_01.scores_0.txt \
             -l 0.001 -e 0.02
    Returns:
        # lambda 0.001                       
        # eta 0.02                           
        4 6 13 17 18 19 20 22 24 26 28 30 49 

    Multi-task without correlation matrix (equivalent to MultiSConES [2]):
    $ python multitask_sfan.py --num_task 2 \
             --networks data/simu1.network.dimacs \
             --nodes_weight data/simu1.scores_0.txt data/simu1.scores_1.txt \
             -l 0.001 -e 0.02 -m 0.01
    $ python multitask_sfan.py --num_tasks 2 \
             --networks ../data/simu_multitask_01.network.dimacs \
             --nodes_weight ../data/simu_multitask_01.scores_0.txt \
                            ../data/simu_multitask_01.scores_1.txt \
             -l 0.001 -e 0.02 -m 0.01 
    Returns:
        # lambda 0.001                           
        # eta 0.02                               
        # mu 0.01                                
        4 6 13 17 18 19 20 22 24 26 28 30 49     
        3 4 7 9 16 19 22 23 27 29 41 43          
    
    Multi-task with correlation matrix:
    $ python multitask_sfan.py --num_tasks 2 \
             --networks ../data/simu_multitask_01.network.dimacs \
             --nodes_weight ../data/simu_multitask_01.scores_0.txt ../data/simu_multitask_01.scores_1.txt \
             --correlation_matrix ../data/simu_multitask_01.task_similarities.txt \
             -l 0.001 -e 0.02 -m 0.01 
    Returns:
        # lambda 0.001
        # eta 0.02
        # mu 0.01
        4 6 13 17 18 19 20 22 24 26 28 30 49
        3 4 7 9 12 16 19 22 23 27 29 41 43  

    References
    ----------
    [1] Azencott, C.-A., Grimm, D., Sugiyama, M., Kawahara, Y., and Borgwardt, K.M. (2013).
    Efficient network-guided multi-locus association mapping with graph cuts.
    Bioinformatics 29, i171--i179.
    [2] Sugiyama, M., Azencott, C., Grimm, D., Kawahara, Y., and Borgwardt, K. (2014).
    Multi-Task Feature Selection on Multiple Networks via Maximum Flows.
    In Proceedings of the 2014 SIAM International Conference on Data Mining, pp. 19--207.
    """
    # Get arguments values
    parser = argparse.ArgumentParser(description = "Generate the super network", add_help = True)
    parser.add_argument("-k", "--num_tasks", help = "Number of tasks", type = int)
    parser.add_argument("-w", "--networks", help = "Paths of networks", nargs = '+')
    parser.add_argument("-r", "--nodes_weight", help = "Paths of node weights", nargs = '+')
    parser.add_argument("-c", "--correlation_matrix", help = "Path of correlation matrix")
    parser.add_argument("-l", "--lambdA", help = "lambda parameter", type = float)
    parser.add_argument("-e", "--eta", help = "eta parameter", type = float)
    parser.add_argument("-m", "--mu", help = "mu parameter", type = float)
    parser.add_argument("-o", "--output", help="File name for runtime output")
    
    args = parser.parse_args()
    num_tasks = args.num_tasks
    
    # Check arguments integrity
    try:
        assert(args.num_tasks >= 1)
    except AssertionError:
        sys.stderr.write("There must be at least one task specified.\n")
        sys.exit(-1)
        
    try:
        assert(len(args.networks) == args.num_tasks or len(args.networks) == 1)
    except AssertionError:
        sys.stderr.write("There must be either 1 network or as many networks as tasks specified.\n")
        sys.exit(-1)
        
    try:
        assert(len(args.nodes_weight) == args.num_tasks)
    except AssertionError:
        sys.stderr.write("There must be as many weight lists as tasks specified.\n")
        sys.exit(-1)
        
    try:
        assert(args.lambdA is not None and args.lambdA > 0.0)
    except AssertionError:
        sys.stderr.write("The lambda parameter must be strictly positive.\n")
        sys.exit(-1)
        
    try:
        assert(args.eta is not None and args.eta > 0.0)
    except AssertionError:
        sys.stderr.write("The eta parameter must be strictly positive.\n")
        sys.exit(-1)
        
    if (args.num_tasks > 1) :
        try:
            assert(args.mu is not None and args.mu > 0.0)
        except AssertionError:
            sys.stderr.write("The mu parameter must be strictly positive.\n")
            sys.exit(-1)

    # Start the process and mark starting time
    time_start = time.clock()

    # Read networks nodes count and edges count
    num_nodes = 0
    num_nodes_each_network = 0
    num_edges = 0
    for task_idx in range(args.num_tasks):
        with open(get_network(args, task_idx), 'r') as f:
            
            ls = f.readline().split()
            
            if not num_nodes_each_network:
                num_nodes_each_network = int(ls[2])
                
            elif (num_nodes_each_network != int(ls[2])) :
                f.close
                sys.stderr.write("All the networks must have the same number of nodes.\n")
                sys.exit(-1)
                
            num_nodes += int(ls[2])
            num_edges += int(ls[3])
            
            f.close()
            
    real_num_nodes = num_nodes + 2
    real_num_edges =  num_edges + (num_nodes_each_network * args.num_tasks) + \
                      (args.num_tasks * (args.num_tasks - 1) * num_nodes_each_network)
            
            
    # Read correlation matrix (if more than one task)
    eta = args.eta
    if (args.num_tasks > 1):
        if args.correlation_matrix:
            # Use the correlation matrix provided
            with open(args.correlation_matrix, 'r') as f:
                correlation_matrix = np.loadtxt(f)
                f.close()
            inv_correlation_matrix = inv(correlation_matrix)
        else:
            # Build the canonical inverse correlation matrix
            eps = args.eta / (10.*args.mu)
            inv_correlation_matrix = -np.ones((args.num_tasks, args.num_tasks)) + \
                                     np.diag(np.ones(args.num_tasks) * (args.num_tasks - 1 + eps) + 1)
            
            # Adjust eta accordingly
            eta = args.eta - args.mu * eps
            try:
                assert(eta > 0)
            except AssertionError:
                sys.stderr.write("Failed to generate a positive eta " + \
                                 "when building the canonical inverse correlation matrix. Check code.\n")
                sys.exit(-1)

        # Sum rows of the correlation matrix
        phi = inv_correlation_matrix.sum(axis=1, dtype='float')

    # Initialize runtimes
    time_task_computations = []
        
    # Generate super-network in dimacs format
    source_node_data = []
    dimacs_graph = ""
    dimacs_graph += ("p max %d %d\n" % (real_num_nodes, real_num_edges))
    dimacs_graph += ("n %d s\n" % (real_num_nodes - 1))
    dimacs_graph += ("n %d t\n" % real_num_nodes)
    
    a = 0.0
    
    # Time count for pre-process step
    time_post_setout_process = time.clock()
    
    for current_task in range(num_tasks):
        f_nt = open(get_network(args, current_task), 'r')
        f_nw = open(args.nodes_weight[current_task], 'r')
        
        ls = f_nt.readline().split()
        while ls[0] != "a":
            ls = f_nt.readline().split()
        
        for i in range(num_nodes_each_network):
            current_node = ((num_nodes_each_network * current_task) + int(ls[1]))
            neighbour_node = ((num_nodes_each_network * current_task) + int(ls[2]))
            
            # Connect nodes within the same task
            while int(ls[1]) == i + 1:
                dimacs_graph += ("a %d %d %f\n" % (current_node, neighbour_node, (float(ls[3])) * args.lambdA))
                ls = f_nt.readline().split()
                if len(ls) == 0:
                    break
                neighbour_node = ((num_nodes_each_network * current_task) + int(ls[2]))
                
            # Connect corresponding nodes across tasks
            if num_tasks > 1:
                target_node = current_node
                while target_node > num_nodes_each_network:
                    target_node -= num_nodes_each_network
                    
                for task_idx in range(num_tasks):
                    real_target_node = num_nodes_each_network * task_idx + target_node
                    
                    if real_target_node != current_node:
                        dimacs_graph += ("a %d %d %f\n" % (current_node, real_target_node, args.mu * \
                                                           inv_correlation_matrix[current_task][task_idx]))
                        
            # Connect nodes to the source and sink nodes
            if num_tasks > 1:
                a = (float(f_nw.readline())) - args.mu * phi[current_task] - eta
            else:
                a = (float(f_nw.readline())) - eta
            
            if a >= 0.0:
                source_node_data.append("a %d %d %f\n" % (real_num_nodes - 1, current_node, a)) 
            else:
                dimacs_graph += ("a %d %d %f\n" % (current_node, real_num_nodes, -a))

        f_nt.close()
        f_nw.close()
        
        time_task_computations.append(time.clock())
        
    for x in source_node_data:
        dimacs_graph += x

    time_all_tasks_computations = time.clock()
        
        
    print "# lambda " + str(args.lambdA)
    print "# eta " + str(args.eta)
    if args.mu:
        print "# mu " + str(args.mu)

    # Pass super network to gt_maxflow
    gt_maxflow.python_entry_point(dimacs_graph, num_nodes_each_network)
    time_total_time = time_gt_maxflow = time.clock()

    sys.stdout.write("\n")
    
    # Process runtimes to a printable string
    runtime_str = ""
    real_time_task_computations = [0 for x in time_task_computations]

    for (i, x) in enumerate(time_task_computations):
        if i > 0:
            real_time_task_computations[i] = x - time_task_computations[i - 1]
            runtime_str += "Task ({0}) computation time: {1}\n".\
                           format(i + 1, x - time_task_computations[i - 1])
            
        else:
            real_time_task_computations[i] = x - time_post_setout_process
            runtime_str += "Task ({0}) computation time: {1}\n".\
                           format(1, x - time_post_setout_process)

    runtime_str += "Task average computation time: {0}\n".\
                   format(np.mean(np.array(real_time_task_computations)))
    runtime_str += "Standard deviation computation time: {0}\n".\
                      format(np.std(np.array(real_time_task_computations)))
    runtime_str += "Network building time: {0}\n".\
                      format(time_all_tasks_computations - time_post_setout_process)
    runtime_str += "gt_maxflow computation time: {0}\n".\
                      format(time_gt_maxflow - time_all_tasks_computations)
    runtime_str += "Process time: {0}\n".format(time_total_time)


    # Save runtime_str to file (if provided), otherwise print to screen
    if args.output is not None:
        with open(args.output, 'w') as output_file:
            output_file.write(runtime_str)
            output_file.close()
    else:
        sys.stdout.write("%s" % runtime_str)
            
        
def get_network(args, task_idx):
    """ Get network file from list of arguments.

    Parameters
    ----------
    args: arguments
        as parsed by argparse

    task_idx: int
        index of the task/network

    Returns
    ------
    Path to the network for task task_idx.    
    """
    if (task_idx == 0 or len(args.networks) == 1) :
        return args.networks[0]
    
    else :
        return args.networks[task_idx]
                
        
if __name__ == "__main__":
    main()

