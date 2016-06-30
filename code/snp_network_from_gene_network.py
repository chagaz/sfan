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

def get_occ_len(occ):
    return occ[2]-occ[1] + 1

def on_same_chromo(occ1, occ2) : 
    return occ1[0] == occ2[0]

def get_overlap(occ1, occ2):
    return max(0, min(occ1[2], occ2[2]) - max(occ1[1], occ2[1]) + 1 )
    

def touch(occ1, occ2):
    if on_same_chromo(occ1, occ2) and 0 < get_overlap(occ1, occ2) < get_occ_len(occ2) : 
        # on same chromo and overlap's length < len(occ2)
        return True
    return False 

def cover(occ1, occ2):
    if on_same_chromo(occ1, occ2) and 0 < get_overlap(occ1, occ2) == get_occ_len(occ2) : 
        # on same chromo and overlap's length == len(occ2)
        return True
    return False 

def get_all_indices(l, e):
    """Get indices of an element in list containing duplicates
    
    Arguments 
    ---------
    l : list 
        can have dupes
    e : an element
        possibly in the list, possibly several times 
    
    Return : 
    --------
    indices : list 
        list of indices where e was seen
    """
    return [index for index in xrange(len(l)) if l[index] == e]
     
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
    genes_occ = dict()
    genes = dict()
    # key : hugo gene symbol 
    # values : tuples containing 3 items : 
    #   - the first = chromo num
    #   - the second = starting position - window value
    #   - the third = ending position + window value
    # genes_occ : one tuple per HUGO gene symbol occurence
    # genes : one tuple per gene duplicates
    # /!\ 'duplicates' =/= occurence of HUGO gene symbol in hugo file
    # see Technotes. TODO : Technotes

    # Save every occurences (exept duplicates) of each Hgs of the file in genes_occ : 
    with open(args.hugo, 'r') as fdHugo:
        # each line : num chromo \t start pos \t end pos \t HUGO gene symbol
        for line in fdHugo:

            line_split = line.split()
            current_Hgs = line_split[3]
            current_data = (int(line_split[0]), int(line_split[1]), int(line_split[2]) ) 

            if current_Hgs not in genes_occ.keys() : 
                genes_occ[current_Hgs] = set() #set not allow duplicate info
 
            genes_occ[current_Hgs].add(current_data)

    
    for Hgs in genes_occ : 

        genes[Hgs] = set()

        # Sort occurences of each Hgs by their length = pos_end - pos_start + 1  
        genes_occ[Hgs] = sorted( genes_occ[Hgs], key = lambda occurence : occurence[2]-occurence[1] + 1 , reverse = True) 

        occ_infos =  [ [None, list()]  for _ in xrange(len( genes_occ[Hgs]) )  ] 
        print "166", occ_infos
        # occ_infos = a list
        # for each occurence : 
        # a tag : 
        #   - None : no info : unseen / seen but nothing to say 
        #   - 'c' : covered
        #   - 't' : touched
        # a list of touched occ
        for current_occ_idx, current_occ in enumerate(genes_occ[Hgs]) : 
            print "175", occ_infos
            if occ_infos[current_occ_idx][0] != 'c' : #if uncovered : 
                # tag info on remaining uncoverded occ:  
                for i, occ in enumerate(genes_occ[Hgs][current_occ_idx+1:]) :
                    occ_idx = i + current_occ_idx + 1
                    print occ_infos
                    if occ_infos[occ_idx][0] != 'c' : 
                        if cover(current_occ, occ): 
                            occ_infos[occ_idx][0] = 'c'
                        elif touch(current_occ, occ) : 
                            occ_infos[occ_idx][0] = 't'
                            occ_infos[current_occ_idx][1].append(occ_idx)

            # otherwise, info is already tagged by previous occ
 

        if (len( genes_occ[Hgs] )) >= 2 : 

            # occ to keep = those that haven't 'c' tag :
            tag_list = [occ_infos[i][0] for i in xrange (len(occ_infos) ) ]
            c_tagged = get_all_indices( tag_list , 'c')

            all_tagged = range(len(occ_infos))
            no_c_tagged = [i for i in reversed(all_tagged) if i not in c_tagged] # a stack-list

            if len (get_all_indices( tag_list , 't') ) > 0 : import pdb; pdb.set_trace() 
            
            while no_c_tagged:
                occ_idx = no_c_tagged.pop()

                current_data = genes_occ[Hgs][occ_idx]

                if not occ_infos[occ_idx][1] : # touches no one
                    genes[Hgs].add(current_data)

                else  : # touches someth
                    touched_list = occ_infos[occ_idx][1]
                    for touched_idx in touched_list : 
                        current_data = get_union(current_data, genes_occ[Hgs][touched_idx])
                        touched_list.extend(occ_infos[touched_idx][1])
                        no_c_tagged.remove(touched_idx)


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
