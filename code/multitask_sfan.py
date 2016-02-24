""" multitask_sfan.py: Generate the super-network for the given multi-task network-guided
feature selection problem and run maxflow on it.

Single task examples, equivalent to SConES [1]:
>>> import subprocess
>>> p = subprocess.Popen(['python', 'multitask_sfan.py', '--num_tasks', '1', \
                           '--networks', '../data/simu_01/simu_01.network.dimacs', \
                           '--node_weights', '../data/simu_01/simu_01.scores_0.txt', \
                           '-l', '0.001', '-e', '0.02', '--output', '/tmp/test'], \
                           stdout=subprocess.PIPE)
>>> p.communicate()[0][:-1]
'# lambda 0.001\\n# eta 0.02\\n4 6 13 17 18 19 20 22 24 26 28 30 49 '

>>> p = subprocess.Popen(['python', 'multitask_sfan.py', '--num_tasks', '1', \
                           '--networks', '../data/simu_01/simu_01.network.dimacs', \
                           '--node_weights', '../data/simu_01/simu_01.scores_1.txt', \
                           '-l', '0.001', '-e', '0.02', '--output', '/tmp/test'], \
                           stdout=subprocess.PIPE)
>>> p.communicate()[0][:-1]
'# lambda 0.001\\n# eta 0.02\\n3 4 7 9 12 16 19 20 22 23 24 27 29 41 43 '

Multi-task with correlation matrix meant to be equivalent to MultiSConES [2]
(Example created with the R code in [2]):
>>> p = subprocess.Popen(['python', 'multitask_sfan.py', '--num_tasks', '2', \
                           '--networks', '../data/simu_ms/simu_ms.network.dimacs', \
                           '--node_weights', '../data/simu_ms/simu_ms.scores_0.txt', \
                                             '../data/simu_ms/simu_ms.scores_1.txt', \
                           '--correlation_matrix', '../data/simu_ms/simu_ms.task_similarities.txt', \
                           '-l', '0.2', '-e', '0.15', '-m', '0.1', \
                           '--output', '/tmp/test'], \
                           stdout=subprocess.PIPE)
>>> p.communicate()[0][:-1]
'# lambda 0.2\\n# eta 0.15\\n# mu 0.1\\n1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 \\n1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 '

Same result without giving the correlation matrix (adjusting the value of eta):
>>> p = subprocess.Popen(['python', 'multitask_sfan.py', '--num_tasks', '2', \
                           '--networks', '../data/simu_ms/simu_ms.network.dimacs', \
                           '--node_weights', '../data/simu_ms/simu_ms.scores_0.txt', \
                                             '../data/simu_ms/simu_ms.scores_1.txt', \
                           '-l', '0.2', '-e', '0.25', '-m', '0.1', \
                           '--output', '/tmp/test'], \
                           stdout=subprocess.PIPE)
>>> p.communicate()[0][:-1]
'# lambda 0.2\\n# eta 0.25\\n# mu 0.1\\n1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 \\n1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 '

The result is different if we do not adjust the value of eta:
>>> p = subprocess.Popen(['python', 'multitask_sfan.py', '--num_tasks', '2', \
                           '--networks', '../data/simu_ms/simu_ms.network.dimacs', \
                           '--node_weights', '../data/simu_ms/simu_ms.scores_0.txt', \
                                             '../data/simu_ms/simu_ms.scores_1.txt', \
                           '-l', '0.2', '-e', '0.15', '-m', '0.1', \
                           '--output', '/tmp/test'], \
                           stdout=subprocess.PIPE)
>>> p.communicate()[0][:-1]
'# lambda 0.2\\n# eta 0.15\\n# mu 0.1\\n1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 1882 1883 1884 1885 1886 1887 1888 1889 1890 1891 1892 \\n1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 '
 
Multitask with correlation matrix:
>>> p = subprocess.Popen(['python', 'multitask_sfan.py', '--num_tasks', '2', \
                           '--networks', '../data/simu_01/simu_01.network.dimacs', \
                           '--node_weights', '../data/simu_01/simu_01.scores_0.txt', \
                                             '../data/simu_01/simu_01.scores_1.txt', \
                           '--correlation_matrix', '../data/simu_01/simu_01.task_similarities.txt', \
                           '-l', '0.001', '-e', '0.02', '-m', '0.01', \
                           '--output', '/tmp/test'], \
                           stdout=subprocess.PIPE)
>>> p.communicate()[0][:-1]
'# lambda 0.001\\n# eta 0.02\\n# mu 0.01\\n4 6 13 17 18 19 20 22 24 26 28 30 49 \\n3 4 7 9 12 16 19 22 23 27 29 41 43 '

Same problem, without correlation matrix:
>>> p = subprocess.Popen(['python', 'multitask_sfan.py', '--num_tasks', '2', \
                           '--networks', '../data/simu_01/simu_01.network.dimacs', \
                           '--node_weights', '../data/simu_01/simu_01.scores_0.txt', \
                                             '../data/simu_01/simu_01.scores_1.txt', \
                           '-l', '0.001', '-e', '0.02', '-m', '0.01', \
                           '--output', '/tmp/test'], \
                           stdout=subprocess.PIPE)
>>> p.communicate()[0][:-1]
'# lambda 0.001\\n# eta 0.02\\n# mu 0.01\\n4 6 13 17 18 19 20 22 24 26 28 30 49 \\n3 4 7 9 12 16 19 20 22 23 24 27 29 41 43 '
    
Same problem, with mu=0 (equivalent to solving both tasks independently, as with SCoNES above):
>>> p = subprocess.Popen(['python', 'multitask_sfan.py', '--num_tasks', '2', \
                           '--networks', '../data/simu_01/simu_01.network.dimacs', \
                           '--node_weights', '../data/simu_01/simu_01.scores_0.txt', \
                                             '../data/simu_01/simu_01.scores_1.txt', \
                           '-l', '0.001', '-e', '0.02', '-m', '0', \
                           '--output', '/tmp/test'], \
                           stdout=subprocess.PIPE)
>>> p.communicate()[0][:-1]
'# lambda 0.001\\n# eta 0.02\\n4 6 13 17 18 19 20 22 24 26 28 30 49 \\n3 4 7 9 12 16 19 20 22 23 24 27 29 41 43 '

References
----------
[1] Azencott, C.-A., Grimm, D., Sugiyama, M., Kawahara, Y., and Borgwardt, K.M. (2013).
    Efficient network-guided multi-locus association mapping with graph cuts.
    Bioinformatics 29, i171--i179.
[2] Sugiyama, M., Azencott, C., Grimm, D., Kawahara, Y., and Borgwardt, K. (2014).
    Multi-Task Feature Selection on Multiple Networks via Maximum Flows.
    In Proceedings of the 2014 SIAM International Conference on Data Mining, pp. 19--207.
"""

import argparse
import doctest
import numpy as np
from numpy.linalg import inv
import sys
import time

import gt_maxflow


def get_network(network_arg, task_idx):
    """ Get network file from list of arguments.

    Parameters
    ----------
    network_arg: list of space-separated network filenames
        as parsed by argparse (--networks)

    task_idx: int
        index of the task/network

    Returns
    ------
    Path to the network for task task_idx.    
    """
    if (task_idx == 0 or len(network_arg) == 1):
        return network_arg[0]
    else:
        return network_arg[task_idx]



class Sfan(object):
    """ Solve a multi-task network-guided feature selection problem.

    Attributes
    ----------
    num_tasks: int
        Number of tasks.
    lbd: float
        Regularization parameter for connectivity.
    eta: float
        Regularization paramter for sparsity.
    networks_f: filename(s)
        List of paths of network files, space-separated.
        If the same network is to be used for all files, only specify one file.
    node_weights_f: filename(s)
        List of paths of network node weights (i.e. feature relevance scores), space-separated.
        One list per task.

    num_nodes_each_network: int
        Number of nodes per task network.
    super_num_nodes: int
        Total number of nodes of the super network.
    super_num_edges: int
        Total number of edges of the super network.

    dimacs_graph: string
        Dimacs description of the super-network.
    
    mu: {float, None}, optional
        Regularization paramter for task relatednes.
    output: {filename, None}, optional
        File where to store computation run times.
    inv_correlation_matrix: {(num_tasks, num_tasks) array, None}, optional
        Inverse of the $\Omega$ matrix of correlations between tasks.
    phi: {(num_tasks, ) array, None}, optional
        Row-wise sum of the $\Omega^{-1}$ inverse matrix of correlations between tasks.

    output_f: {filename, None}, optional
        File where to store computation run times.        
    """
    def __init__(self, num_tasks, networks_f, node_weights_f, lbd, eta, mu=None,
                 correlation_matrix_f=None, output_f=None):
        """
        Parameters
        ----------
        num_tasks: int
            Number of tasks.
        networks_f: filename(s)
            List of paths of network files, space-separated.
            If the same network is to be used for all files, only specify one file.
        node_weights_f: filename(s)
            List of paths of network node weights (i.e. feature relevance scores), space-separated.
            One list per task.
        lbd: float
            Regularization parameter for connectivity.
        eta: float
            Regularization paramter for sparsity.
        mu: {float, None}, optional
            Regularization paramter for task relatednes.
        correlation_matrix_f: {filename, None}, optional
            Path of correlation matrix.
        output_f: {filename, None}, optional
            File where to store computation run times.        
        """
        self.num_tasks = num_tasks
        self.networks_f = networks_f
        self.node_weights_f = node_weights_f
        self.lbd = lbd
        self.eta = eta
        self.mu = mu
        self.output_f = output_f
        
        # Read networks nodes count and edges count
        self.super_num_nodes = 0
        self.num_nodes_each_network = 0
        self.super_num_edges = 0
        for task_idx in range(num_tasks):
            with open(get_network(networks_f, task_idx), 'r') as f:
                ls = f.readline().split()

                if not self.num_nodes_each_network:
                    self.num_nodes_each_network = int(ls[2])

                elif (self.num_nodes_each_network != int(ls[2])) :
                    f.close
                    sys.stderr.write("All the networks must have the same number of nodes.\n")
                    sys.exit(-1)

                self.super_num_nodes += int(ls[2])
                self.super_num_edges += int(ls[3])

                f.close()

        # The super network has one node for each node in the task networks, + source and sink
        self.super_num_nodes += 2
        # Edges in the super-network are either:
        #    - initial connections inside each network
        #    - connections of each node to source or sink (directed)
        #    - connections between the corresponding nodes of each network (undirected)
        self.super_num_edges =  self.super_num_edges + \
                                (self.num_nodes_each_network * self.num_tasks) + \
                                (self.num_tasks * (self.num_tasks - 1) * self.num_nodes_each_network)
        
        # Read correlation matrix (if more than one task)
        if (self.num_tasks > 1):
            if correlation_matrix_f:
                # Use the correlation matrix provided
                with open(correlation_matrix_f, 'r') as f:
                    correlation_matrix = np.loadtxt(f)
                    f.close()
                self.inv_correlation_matrix = inv(correlation_matrix)
            else:
                # Build the canonical inverse correlation matrix
                if self.mu > 0:
                    eps = self.eta / (10.*self.mu)
                else:
                    eps = 0.01 * self.eta
                self.inv_correlation_matrix = -np.ones((self.num_tasks, self.num_tasks)) + \
                                              np.diag(np.ones(self.num_tasks) * \
                                                      (self.num_tasks - 1 + eps) + 1)

                # Adjust eta accordingly
                self.eta = self.eta - self.mu * eps
                try:
                    assert(self.eta > 0)
                except AssertionError:
                    sys.stderr.write("Failed to generate a positive eta " + \
                                     "when building the canonical inverse correlation matrix. Check code.\n")
                    sys.exit(-1)

            # Sum rows of the correlation matrix
            self.phi = self.inv_correlation_matrix.sum(axis=1, dtype='float')

            # print self.inv_correlation_matrix
            # sys.exit(0)

            
    def create_dimacs(self):
        """ Create the dimacs description of the problem's super-network.

        Modified attributes
        -------------------
        self.dimacs_graph: string
            Dimacs description of the super-network.

        Return
        ------
        time_task_computations: list
            Computation times for each task.
        """
        # Initialize runtimes
        time_task_computations = []
        
        source_node_data = []
        self.dimacs_graph = ""
        self.dimacs_graph += ("p max %d %d\n" % (self.super_num_nodes, self.super_num_edges))
        self.dimacs_graph += ("n %d s\n" % (self.super_num_nodes - 1))
        self.dimacs_graph += ("n %d t\n" % self.super_num_nodes)

        a = 0.0

        for current_task in range(self.num_tasks):
            with open(get_network(self.networks_f, current_task), 'r') as f_nt:
                with open(self.node_weights_f[current_task], 'r') as f_nw:
                    ls = f_nt.readline().split()
                    while ls[0] != "a":
                        ls = f_nt.readline().split()

                    for node_idx in range(self.num_nodes_each_network):
                        current_node = ((self.num_nodes_each_network * current_task) + int(ls[1]))
                        neighbour_node = ((self.num_nodes_each_network * current_task) + int(ls[2]))

                        # Connect nodes within the same task
                        while int(ls[1]) == node_idx + 1:
                            self.dimacs_graph += ("a %d %d %f\n" % (current_node, neighbour_node,
                                                               (float(ls[3])) * self.lbd))
                            ls = f_nt.readline().split()
                            if len(ls) == 0:
                                break
                            neighbour_node = ((self.num_nodes_each_network * current_task) + int(ls[2]))

                        # Connect corresponding nodes across tasks
                        if self.num_tasks > 1:
                            target_node = current_node
                            while target_node > self.num_nodes_each_network:
                                target_node -= self.num_nodes_each_network

                            for task_idx in range(self.num_tasks):
                                real_target_node = self.num_nodes_each_network * task_idx + target_node

                                if real_target_node != current_node:
                                    self.dimacs_graph += ("a %d %d %f\n" % \
                                                     (current_node, real_target_node,
                                                      self.mu * \
                                                      self.inv_correlation_matrix[current_task][task_idx]))

                        # Connect nodes to the source and sink nodes
                        if self.num_tasks > 1:
                            a = (float(f_nw.readline())) - self.mu * self.phi[current_task] - self.eta
                        else:
                            a = (float(f_nw.readline())) - self.eta

                        if a >= 0.0:
                            source_node_data.append("a %d %d %f\n" % \
                                                    (self.super_num_nodes-1, current_node, a)) 
                        else:
                            self.dimacs_graph += ("a %d %d %f\n" % (current_node, self.super_num_nodes, -a))
                    f_nw.close()
                f_nt.close()

            time_task_computations.append(time.clock())

        for x in source_node_data:
            self.dimacs_graph += x

        return time_task_computations

        
    def run_maxflow(self):
        """ Run gt_maxflow on the super-network.

        Border effect
        -------------
        Prints to screen the list of nodes selected in each network. 
        """
        # Pass super network to gt_maxflow
        gt_maxflow.python_entry_point(self.dimacs_graph, self.num_nodes_each_network)


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
    args.node_weight: filename(s)
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
    Also see: test_multitask_sfan.sh

    $ python multitask_sfan.py --num_tasks 2 \
             --networks ../data/simu_01/simu_01.network.dimacs \
             --node_weights ../data/simu_01/simu_01.scores_0.txt ../data/simu_01/simu_01.scores_1.txt \
             --correlation_matrix ../data/simu_01/simu_01.task_similarities.txt \
             -l 0.001 -e 0.02 -m 0.01
    Returns:
        # lambda 0.001
        # eta 0.02
        # mu 0.01
        4 6 13 17 18 19 20 22 24 26 28 30 49
        3 4 7 9 12 16 19 22 23 27 29 41 43
    """
    # Get arguments values
    parser = argparse.ArgumentParser(description="Generate the super network", add_help=True)
    parser.add_argument("-k", "--num_tasks", help="Number of tasks", type=int)
    parser.add_argument("-w", "--networks", help="Paths of networks", nargs='+')
    parser.add_argument("-r", "--node_weights", help="Paths of node weights", nargs='+')
    parser.add_argument("-t", "--test", help="Run tests", action='store_true')
    parser.add_argument("-c", "--correlation_matrix", help="Path of correlation matrix")
    parser.add_argument("-l", "--lambdA", help="lambda parameter", type=float)
    parser.add_argument("-e", "--eta", help="eta parameter", type=float)
    parser.add_argument("-m", "--mu", help="mu parameter", type=float)
    parser.add_argument("-o", "--output", help="File name for runtime output")
    
    args = parser.parse_args()
    num_tasks = args.num_tasks

    # Testing (with docstring)
    if args.test:
        sys.stdout.write("If no output, all tests passed!\n")
        doctest.testmod()
        sys.exit(0)

    # Check arguments integrity
    try:
        assert(args.num_tasks >= 1)
    except AssertionError:
        sys.stderr.write("There must be at least one task specified.\n")
        sys.stderr.write("Use --help for help.\n")
        sys.exit(-1)
        
    try:
        assert(len(args.networks) == args.num_tasks or len(args.networks) == 1)
    except AssertionError:
        sys.stderr.write("There must be either 1 network or as many networks as tasks specified.\n")
        sys.stderr.write("Use --help for help.\n")
        sys.exit(-1)
        
    try:
        assert(len(args.node_weights) == args.num_tasks)
    except AssertionError:
        sys.stderr.write("There must be as many weight lists as tasks specified.\n")
        sys.stderr.write("Use --help for help.\n")
        sys.exit(-1)
        
    try:
        assert(args.lambdA is not None and args.lambdA > 0.0)
    except AssertionError:
        sys.stderr.write("The lambda parameter must be strictly positive.\n")
        sys.stderr.write("Use --help for help.\n")
        sys.exit(-1)
        
    try:
        assert(args.eta is not None and args.eta > 0.0)
    except AssertionError:
        sys.stderr.write("The eta parameter must be strictly positive.\n")
        sys.stderr.write("Use --help for help.\n")
        sys.exit(-1)
        
    if (args.num_tasks > 1) :
        try:
            assert(args.mu is not None and args.mu >= 0.0)
        except AssertionError:
            sys.stderr.write("The mu parameter must be strictly positive.\n")
            sys.stderr.write("Use --help for help.\n")
            sys.exit(-1)

    # Time stamp: beginning of computations
    time_start = time.clock()

    # Instantiate a sfan solver
    sfan_solver = Sfan(args.num_tasks, args.networks, args.node_weights, args.lambdA, args.eta,
                       mu=args.mu, correlation_matrix_f=args.correlation_matrix, output_f=args.output)

    # Time stamp: end of preprocessing
    time_post_setout_process = time.clock()
        
    # Generate super-network in dimacs format
    time_task_computations = sfan_solver.create_dimacs()

    # Time stamp: end of generation of the super-network
    time_all_tasks_computations = time.clock()
        
    print "# lambda " + str(args.lambdA)
    print "# eta " + str(args.eta)
    if args.mu:
        print "# mu " + str(args.mu)

    # Solve optimization problem (run maxflow)
    sfan_solver.run_maxflow()
    sys.stdout.write("\n")

    # Time stamp: end of optimization
    time_total_time = time_gt_maxflow = time.clock()
    
    # Process runtimes into a printable string
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
    if sfan_solver.output_f is not None:
        with open(sfan_solver.output_f, 'w') as output_file:
            output_file.write(runtime_str)
            output_file.close()
    else:
        sys.stdout.write("%s" % runtime_str)
                
        
if __name__ == "__main__":
    main()

