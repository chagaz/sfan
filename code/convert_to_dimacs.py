# Convert W.txt to dimacs format (adding source and sink)
#
# chloe-agathe.azencott@mines-paristech.fr

import argparse
import numpy as np
import os
import subprocess
import sys

cur_dir = os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0])))

def main():
    """
    Create .dimacs file from W.txt
    """
    parser = argparse.ArgumentParser(description='Create dimacs file', add_help=True)
    parser.add_argument("dir_name", help="name of directory containing the simulated data")
    args = parser.parse_args()

    data_dir = '%s/data/%s' % (cur_dir, args.dir_name)   

    # number of SNPs
    readme_f = "%s/readme.txt" % data_dir
    with open(readme_f, 'r') as f:
        for line in f:
            ls = line.split()
            if ls[1] == 'features':
                num_features = int(ls[0])
        f.close()

    # # random source to node or node to sink capacities
    # capacities = 2*np.random.random(num_features)-1
        
    Wf = "%s/W.txt" % data_dir
    dimacs_f = "%s/W.dimacs" % data_dir

    # number of edges
    cmd = "wc -l %s" % Wf
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    (stdo, stde) = proc.communicate()
    num_edges = int(stdo.split()[0]) #+ num_features
    
    with open(dimacs_f, 'w') as g:
        g.write("p max %d %d\n" % ((num_features+2), num_edges))
        g.write("n %d s\n" % (num_features+1))
        g.write("n %d t\n" % (num_features+2))
        # for (idx, edge) in enumerate(capacities):
        #     if edge > 0:
        #         g.write("a %d %d %.2f\n" % ((num_features+1), (idx+1), edge))
        #     else:
        #         g.write("a %d %d %.2f\n" % ((idx+1), (num_features+2), -edge))
        with open(Wf, 'r') as f:
            for line in f:
                ls = line.split()
                u = int(ls[0])
                v = int(ls[1])
                if u < v:                
                    #g.write("a %s" % line)
                    if len(ls) > 2:
                        g.write("a %d %d %s\n" % (u, v, ls[2]))
                        g.write("a %d %d %s\n" % (v, u, ls[2]))
                    else:
                        g.write("a %d %d 1\n" % ((u+1), (v+1)))
                        g.write("a %d %d 1\n" % ((v+1), (u+1)))                        
            f.close()
        g.close()

    sys.stdout.write("Wrote %s\n" % dimacs_f)
            

if __name__ == '__main__':
    main()
