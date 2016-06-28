# snp_network_from_gene_network.py -- Create SNPs network
#
# jean-daniel.granet@mines-paristech.fr

import sys
import argparse
import time
from operator import itemgetter
import numpy as np # numerical python module
import scipy.sparse as sp # scientific python sparse module

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
    # read hugogenes.txt and save into a dictionnary : the key is the name of the gene
    print 'Save the genes positions : ',
    Start = time.time()
    genes = dict()
    with open(args.hugo, 'r') as fdHugo:
        # each line : num chromo \t start pos \t end pos \t HUGO gene symbol
        for line in fdHugo:
            line_split = line.split()
            if line_split[3] in genes.keys(): # handling HUGO gene symbols duplicates
                if genes[line_split[3]][1] > int(line_split[1]) - int(args.window):
                    genes[line_split[3]][1] = int(line_split[1]) - int(args.window)
                elif genes[line_split[3]][2] > int(line_split[2]) + int(args.window):
                    genes[line_split[3]][2] = int(line_split[2]) - int(args.window)
            else:
                genes[line_split[3]] = [int(line_split[0]), int(line_split[1]) - int(args.window), int(line_split[2]) + int(args.window)]
        fdHugo.close()
    print '\033[92m' + 'DONE' + '\033[0m'
    End = time.time()
    print 'Exec time :' + str(End - Start)
    #---------------------------------------------------------------------------
    
    #---------------------------------------------------------------------------
    # attach the SNPs to the genes
    print 'Attach each SNPs of a gene to each other : ',
    Start = time.time()
    for keyG in genes:
        for idxSNP in xrange(Chromosome[genes[keyG][0] - 1][0], Chromosome[genes[keyG][0] - 1][1] + 1):
            if SNPs[idxSNP][0] == genes[keyG][0]:
                if SNPs[idxSNP][1] > genes[keyG][1] and SNPs[idxSNP][1] < genes[keyG][2]:
                    genes[keyG].append(idxSNP)
        for idxG in xrange(3, len(genes[keyG])):
            for idxG2 in xrange(idxG + 1, len(genes[keyG])):
                if net[genes[keyG][idxG2], genes[keyG][idxG]] != 1:
                    net[genes[keyG][idxG], genes[keyG][idxG2]] = 1
    print '\033[92m' + 'DONE' + '\033[0m'
    End = time.time()
    print 'Exec time :' + str(End - Start)
    #---------------------------------------------------------------------------
    
    #---------------------------------------------------------------------------
    # connect the SNPs of gene A to the SNPs of gene B
    print 'Connect each SNPs of gene A to each SNPs of gene B : ',
    Start = time.time()
    with open(args.acsn, 'r') as fdAcsn:
        for line in fdAcsn:
            line_split = line.split()
            if line_split[0] in genes and line_split[2] in genes:
                if len(genes[line_split[0]]) > 3 and len(genes[line_split[2]]) > 3:
                    for idxGA in xrange(3, len(genes[line_split[0]])):
                        for idxGB in xrange(3, len(genes[line_split[2]])):
                            if net[genes[line_split[2]][idxGB], genes[line_split[0]][idxGA]] != 1:
                                net[genes[line_split[0]][idxGA], genes[line_split[2]][idxGB]] = 1
        fdAcsn.close()
    print '\033[92m' + 'DONE' + '\033[0m'
    End = time.time()
    print 'Exec time :' + str(End - Start)
    #---------------------------------------------------------------------------
    
    #---------------------------------------------------------------------------
    # write the network into the output file
    print 'Write the network into the output file : ',
    Start = time.time()
    [X, Y] = net.nonzero()
    array_xy = zip(X, Y)
    with open(args.output, 'w') as fdOutput:
        for (x, y) in array_xy:
            if x != y:
                fdOutput.write('%s\n' % '\n'.join(['%s %s %s %s' % (SNPs[x][0], SNPs[x][1], SNPs[y][0], SNPs[y][1])]))
        fdOutput.close()
    print '\033[92m' + 'DONE' + '\033[0m'
    End = time.time()
    print 'Exec time :' + str(End - Start)
    #---------------------------------------------------------------------------
    
    
    TotalEndTime = time.time()
    print "Total excution time :" + str(TotalEndTime - TotalStartTime)

if __name__ == "__main__":
    main()
