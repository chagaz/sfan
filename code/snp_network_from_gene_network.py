# snp_network_from_gene_network.py -- Create SNPs network
#
# jean-daniel.granet@mines-paristech.fr

import sys
import argparse
import time
from operator import itemgetter
import numpy as np # numerical python module
import scipy.sparse as sp # scientific python sparse module
from sympy import Interval, Union

def main():
    
    TotalStartTime = time.time()
    
    #---------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Create SNPs network')
    parser.add_argument('acsn', help='gene network')
    parser.add_argument('map', help='SNPs positions')
    parser.add_argument('hugo', help='gene positions')
    parser.add_argument('window', help='window for SNP-gene association')
    parser.add_argument('output', help='output file')
    args = parser.parse_args()
    #---------------------------------------------------------------------------
    
    #---------------------------------------------------------------------------
    # create a list wich contains SNPs positions. Using PLINK's internal numeric coding for chromosome number
    print 'Creation of the list of SNPS positions : ',
    Start = time.time()
    SNPs = list()
    with open(args.map, 'r') as fdMap:
        # Map file structure description : http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map
        for (line_idx, line) in enumerate(fdMap):
            if line_idx > 0: # avoid header of map file
                line_split = line.split()
                if int(line_split[3]) >= 0: # "To exclude a SNP from analysis, set the 4th column (physical base-pair position) to any negative value" -- http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map
                    if line_split[0] == 'X': # X chromosome
                        SNPs.append((23, int(line_split[3])))
                    elif line_split[0] == 'Y': # Y chromosome
                        SNPs.append((24, int(line_split[3])))
                    elif line_split[0] == 'XY': # Pseudo-autosomal region of X
                        SNPs.append((25, int(line_split[3])))
                    elif line_split[0] == 'MT': # Mitochondrial
                        SNPs.append((26, int(line_split[3])))
                    else:
                        SNPs.append((int(line_split[0]), int(line_split[3])))
        fdMap.close()
    End = time.time()
    print '\033[92m' + 'DONE' + '\033[0m'
    print 'Exec time :' + str(End - Start)
    #---------------------------------------------------------------------------
    
    # sort the SNPs array by chromosome
    SNPs.sort(key=itemgetter(0, 1))
    
    # create a sparse matrix of size len(SNPs) x len(SNPs) to create the network
    print 'Creation of the matrix to save the network : ',
    net = sp.lil_matrix((len(SNPs), len(SNPs)))
    print '\033[92m' + 'DONE' + '\033[0m'

    #---------------------------------------------------------------------------
    # connect each SNPs to the nearest
    print 'Connect each SNPs to the nearest : ',
    Start = time.time()
    
    Chromosome = [(-1,-1)] * 26 # 26 = 22 pairs of autosomes (1-22) + 3 gonosomes : X (23), Y (24), XY (25) + 1 mitochondrial : MT (26)
    # This list saves, for each chromo, which SNP is at the begining, which one is at the end
    
    for idxSNP in xrange(0, len(SNPs)):
        # handle first read SNPs that can not be on the first Chromo
        if idxSNP == 0 : 
            Chromosome[SNPs[idxSNP][0] -1 ] = (0, -1)  #/!\ -1 needed because in Python, chromo numeratation begin to 0, not 1
            # Can't do Chromosome[SNPs[idxSNP][0][0] because tuple do not support item assignment
        if idxSNP + 1 < len(SNPs): # the current SNP is not the last. +1 is needed because numerotation begin to 0
            if SNPs[idxSNP][0] != SNPs[idxSNP + 1][0]: # the current SNP and the next one are not on the same chromosome
                # save the current SNP_index as the last SNP of the current chromo : 
                Chromosome[SNPs[idxSNP][0] - 1] = (Chromosome[SNPs[idxSNP][0] - 1][0], idxSNP) # /!\ -1 needed because in Python, chromo numeratation begin to 0, not 1
                # save the next SNP index as the first of the next chromo : 
                Chromosome[SNPs[idxSNP][0]] = (idxSNP + 1, 0)
            if SNPs[idxSNP][0] == SNPs[idxSNP + 1][0]: # the current SNP and the next one are on the same chromosome -> connect them
                net[idxSNP, idxSNP + 1] = 1
        else: # the current studied SNP is the last one
            Chromosome[SNPs[idxSNP][0] - 1] = (Chromosome[SNPs[idxSNP][0] - 1][0], idxSNP)
    print '\033[92m' + 'DONE' + '\033[0m'
    End = time.time()
    print 'Exec time :' + str(End - Start)
    #---------------------------------------------------------------------------

    
    #---------------------------------------------------------------------------
    # read hugogenes.txt and save gene duplicates of each Hgs of the file in
    # into a dictionnary : 
    print 'Save the genes positions : ',
    Start = time.time()
    genes = dict()
    # key : hugo gene symbol 
    # values : tuples containing 2 items : 
    #   - the first = chromo num
    #   - the second = Interval 
    #     ( starting position - window value, ending position + window value)
    #   - the third = list of SPNs (empty list at this steps)
    # genes : one tuple per chromosome where duplicates are founded
    # /!\ 'genes duplicates' =/= occurence of HUGO gene symbol in hugo file


    with open(args.hugo, 'r') as fdHugo:
        # each line : num chromo \t start pos \t end pos \t HUGO gene symbol
        for line in fdHugo:

            # get data from file : 
            line_split = line.split()
            current_Hgs = line_split[3]
            current_chromo_num = int(line_split[0])
            current_Interval = Interval(int(line_split[1]), int(line_split[2]) )
            current_data = (current_chromo_num, current_Interval , list() ) 


            if current_Hgs not in genes.keys() : 
                genes[current_Hgs] = list() 

            print 'before', genes[current_Hgs]
            print current_data

            # handle multi occurence and save new data
            chromo_num_list = [genes[current_Hgs][i][0] for i in xrange (len (genes[current_Hgs]))]
            if current_chromo_num in chromo_num_list : # there is possibly an overlap 
                # get the index of the tuple holding info on the same chromo...
                duplicate_idx = chromo_num_list.index(current_chromo_num) 
                # ... and the associated data, 
                duplicate_data = genes[current_Hgs][duplicate_idx]
                # merge old and current data, 
                to_save = (current_chromo_num, Union (current_Interval, duplicate_data[1]), list() )
                # and save the merged data : 
                genes[current_Hgs][duplicate_idx] = to_save
            else : # no overlap are possible,
                # thus just add a new tuple holding current data : 
                genes[current_Hgs].append(current_data)

            print 'after', genes[current_Hgs]
            print '-------------------'

    print '\033[92m' + 'DONE' + '\033[0m'
    End = time.time()
    print 'Exec time :' + str(End - Start)
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    # attach the SNPs to the genes

    print 'Attach each SNPs of a gene to each other : ',
    Start = time.time()
    for Hgs in genes:

        #----------
        # List SNPs belonging to a gene :
        SNPs_in_Hgs = list() 
        # Each Hgs can be on several chromo. 
        # Get the list of chromo num of current Hgs: 
        list_chromo_num = [dupe[0]-1 for dupe in genes[Hgs] ]
        # /!\ -1 because numeration begin to 0 and not to 1

        # For each chromo num, we have an union of interval : 

        for dupe_idx, chromo_num in enumerate(list_chromo_num) : 
            # Get the range of SNP that are positionated on this chromo : 
            num_SNP_range = xrange(Chromosome[chromo_num][0], Chromosome[chromo_num][1] + 1)
            # For each SNP in this range : 
            for SNP_idx in num_SNP_range : 
                # If the SNP belong to gene dupe on this chromo (take the window into account)
                if SNPs[SNP_idx][1] in genes[Hgs][dupe_idx][1] : 
                    # Add the gene to the list : 
                    SNPs_in_Hgs.append(SNP_idx)
                    genes[Hgs][dupe_idx][2].append(SNP_idx)
        #----------
        # Attach each SNPs of a gene to each other :
        for SNP_idx1 in xrange(len(SNPs_in_Hgs)):
            for SNP_idx2 in xrange(SNP_idx1 + 1, len(SNPs_in_Hgs)):
                if net[SNPs_in_Hgs[SNP_idx2], SNPs_in_Hgs[SNP_idx2]] != 1: #why ??? 
                    net[SNPs_in_Hgs[SNP_idx1], SNPs_in_Hgs[SNP_idx2]] = 1
                else : import pdb; pdb.set_trace() 
        print 'bim'
    print '\033[92m' + 'DONE' + '\033[0m'
    End = time.time()
    print 'Exec time :' + str(End - Start)
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    # connect the SNPs of gene A to the SNPs of gene B
    print 'Connect each SNP of Hgs A to each SNP of Hgs B : ',
    Start = time.time()
    with open(args.acsn, 'r') as fdAcsn:
        # each line : <hgsA> <name of relationship> <hgsB>
        for line in fdAcsn:
            line_split = line.split()
            HgsA = line_split[0]
            HgsB = line_split[2] 
            # we can only use Hgs for which we have infos about localisation : 
            if HgsA in genes and HgsB in genes:
                # get SNPs of these Hgs : 
                SNPs_of_A = [SNP_idx for dupe_infos in genes[HgsA] for SNP_idx in dupe_infos[2]  ]
                SNPs_of_B = [SNP_idx for dupe_infos in genes[HgsB] for SNP_idx in dupe_infos[2]  ]
                # if these 2 interacting Hgs have some SNPs : 
                # genes[Hgs] = (chromo num, Interval, [list, of, SNPs, indices])
                if len( SNPs_of_A ) > 0 and len(SNPs_of_B) > 0:
                    # Connect each SNP of hgsA to each SNP of hgsB : 
                    for SNPA in SNPs_of_A:
                        for SNPB in SNPs_of_B:
                            if net[SNPB, SNPA] != 1: # why ???
                                net[ SNPA, SNPB] = 1
                            else : import pdb; pdb.set_trace()
        fdAcsn.close()
    print '\033[92m' + 'DONE' + '\033[0m'
    End = time.time()
    print 'Exec time :' + str(End - Start)
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    # write the network into the output file
    print 'Write the network into the output file : ',
    # dimacs format

    Start = time.time()

    [X, Y] = net.nonzero()
    #X : line number with non zero data
    #Y : for each X, col number with non zero data
    # /!\ dimacs numeration begin to 1 but python to 0 !! 
    # so move the numeratation : 
    X += 1 
    Y += 1
    array_xy = zip(X, Y)
    # array_xy is a list of tuples, where i-th tuple contains the i-th element from X and Y

    with open(args.output, 'w') as fdOutput:
        # write problem line : 
        fdOutput.write('p max %d %d \n' % (len(SNPs), len(array_xy)*2 ) ) # * 2 because edges are written in both direction
        # write an Arc Descriptor
        # for each coord of non zero value in net : 
        for (x, y) in array_xy:
            if x != y: # if non zero value is not on diag... why ??? 
                # write 
                fdOutput.write('a %d %d 1\n' % (x, y) )
                fdOutput.write('a %d %d 1\n' % (y, x) )
                # Remark : need an undirected graph
                # -> we give the link in both direction
        fdOutput.close()
    print '\033[92m' + 'DONE' + '\033[0m'
    End = time.time()
    print 'Exec time :' + str(End - Start)
    #---------------------------------------------------------------------------
    
    
    TotalEndTime = time.time()
    print "Total excution time :" + str(TotalEndTime - TotalStartTime)

if __name__ == "__main__":
    main()
