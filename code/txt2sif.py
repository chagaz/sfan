# txt2sif.py -- Convert .txt into a .sif format


import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Convert .txt to .sif')
    parser.add_argument('txt', help="""txt file to convert. Fields are : 
<gene 1 name in ACSN> <name of relationship> <gene 2 name in ACSN> <semi-colon separated list of PubMedID >. Columns are TAB separated
    """)
    parser.add_argument('sif', help='SIF file to create.')
    args = parser.parse_args()
    
    # read the txt file and convert it into a sif file
    with open(args.txt, 'r') as fdr:
        with open(args.sif, 'w') as fdw:
            for line in fdr:
                line_split = line.split('\t')
                fdw.write("%s\t%s\t%s\n" % (line_split[0], line_split[1], line_split[2]))
            fdw.close()
        fdr.close()

if __name__ == "__main__":
    main()
