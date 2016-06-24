# csv2map.py -- Convert .csv into a .map format
#   Description of MAP format: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map
#
# jean-daniel.granet@mines-paristech.fr

import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Convert .csv to .map')
    parser.add_argument('csv', help="""CSV file to convert. Fields are : Index,Illumina_SNP_Name,Alternative_SNP_Name,Chromosome,Build36_Position,Build37_Position,new_rsname,Strand,TopAlleles,ForwardAlleles,DesignAlleles
    """)
    parser.add_argument('map', help='MAP file to create')
    args = parser.parse_args()
    
    # read the csv file and convert it into a MAP file
    with open(args.csv, 'r') as fdr:
        with open(args.map, 'w') as fdw:
            for line in fdr:
                line_split = line.split(',')
                fdw.write("%s\n" % "\n".join(["%s %s 0 %s" % (line_split[3], line_split[2], line_split[5])]))
            fdw.close()
        fdr.close()

if __name__ == "__main__":
    main()
