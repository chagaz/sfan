# acsn2hugo.py -- Convert ACSN gene to HUGO gene
#
# jean-daniel.granet@mines-paristech.fr

import sys
import argparse

def main():
    # get arguments : 
    parser = argparse.ArgumentParser(description='Convert ACSN gene to HUGO gene')
    parser.add_argument('acsn', help='ACSN to convert in .sif format')
    parser.add_argument('curated', help='reference file in gmt format')
    parser.add_argument('hugo', help='HUGO file to create')
    args = parser.parse_args()
    
    curated = dict()
    # read the curated file and save into a dictionary :
    # the key is the acsn gene name and the value is the hugo gene name
    # curated [ "acsn entity name" ] = "hugo\tgene\tnames\tin\tthe\tacsn\tentity"
    with open(args.curated, 'r') as fdCurated:
        for line in fdCurated:
            line_split = line.split('\tna\t')
            curated[line_split[0].strip()] = line_split[1].strip()
        fdCurated.close()

    # read the acsn network file and write the corresponding hugo network into the output file
    with open(args.acsn, 'r') as fdAcsn:
        with open(args.hugo, 'w') as fdHugo:
            for line in fdAcsn:
                line_split = line.split('\t')
                geneA = curated[line_split[0].strip()].split()
                geneB = curated[line_split[2].strip()].split()
                for elemA in geneA:
                    for elemB in geneB:
                        fdHugo.write("%s\n" % "\n".join(["%s %s %s" % (elemA, line_split[1], elemB)]))
            fdHugo.close()
        fdAcsn.close() 

if __name__ == "__main__":
    main()
