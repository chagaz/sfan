# acsn2hugo.py -- Convert ACSN entities network to HUGO genes network
#
# jean-daniel.granet@mines-paristech.fr

import sys
import argparse

def main():
    # get arguments : 
    parser = argparse.ArgumentParser(description='Convert ACSN entities network to HUGO genes network')
    parser.add_argument('acsn_network_fname', help='ACSN entities network to convert, in .sif format')
    parser.add_argument('correspondances_fname', help='reference file : correspondance between ACSN entities and HUGO genes symbols, in .gmt format')
    parser.add_argument('output_fname', help='file to create : genes network using HUGO gene symbols, in sif format')
    args = parser.parse_args()
    
    correspondances = dict()
    # read the correspondances file and save into a dictionary :
    # the key is the acsn gene name and the value is the hugo gene name
    # correspondances [ "acsn entity name" ] = "hugo\tgene\tnames\tin\tthe\tacsn\tentity"
    with open(args.correspondances_fname, 'r') as fdCorrespondances:
        # File structure : acsn entity name\t description (=na)\tgene\list\in\tHUGO\symbol
        for line in fdCorrespondances:
            line_split = line.split('\tna\t')
            ACSN_entity = line_split[0].strip()
            HUGO_gene_symbols_in_ACSN_entity = line_split[1].strip().split('\t')
            
            if ACSN_entity not in correspondances : 
                correspondances[ACSN_entity] = set() 
            
            correspondances[ACSN_entity].update(HUGO_gene_symbols_in_ACSN_entity) #do not keep duplicates
        fdCorrespondances.close()

    # read the acsn network file and write the corresponding hugo network into the output file
    with open(args.acsn_network_fname, 'r') as fdAcsn:
        # File structure : acsn entity name 1\tName of relationship\t acsn entity name 2 (\t PMIDS;PMIDS;PMIDS...)
        # even if not sif and has header : no problem since no corresponding HUGO symbol will be found
        with open(args.output_fname, 'w') as fdOut:
            #File structure : hugo gene symbol\tName of relationship\thugo gene symbol
            for line in fdAcsn: # "acsn entity name 1\tName of relationship\t acsn entity name 2"
                line_split = line.split('\t')# ["acsn entity name 1","Name of relationship","acsn entity name 2"]
                ACSN_entity_A = line_split[0].strip()
                ACSN_entity_B = line_split[2].strip()
                HUGO_genes_symbols_in_ACSN_entity_A = correspondances.get(ACSN_entity_A)
                HUGO_genes_symbols_in_ACSN_entity_B = correspondances.get(ACSN_entity_B)
                if HUGO_genes_symbols_in_ACSN_entity_A and HUGO_genes_symbols_in_ACSN_entity_B : 
                    for elemA in HUGO_genes_symbols_in_ACSN_entity_A:
                        for elemB in HUGO_genes_symbols_in_ACSN_entity_B:
                            fdOut.write("%s\n" % "\n".join(["%s %s %s" % (elemA, line_split[1], elemB)]))
            fdOut.close()
        fdAcsn.close() 

if __name__ == "__main__":
    main()
