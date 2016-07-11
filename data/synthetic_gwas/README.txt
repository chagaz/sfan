FILE FORMAT

acsn_names.gmt:
    Correspondence between ACSN entities and HUGO names (gmt format)
    From: https://acsn.curie.fr/files/acsn_names.gmt
	Each line corresponds to one ACSN entity
	<entity name in ACSN> <dataset name> <gene name 1 in HUGO>  <gene name 2 in HUGO> ...
	columns are TAB separated
	each column can contain SPACES

acsn_ppi_ver2.txt:
    PPI interactions (Tab-delimited text format)
    From: https://acsn.curie.fr/files/acsn_ppi_ver2.txt
	Each line corresponds to one interaction (network edge)
	<gene 1 name in ACSN> <name of relationship> <gene 2 name in ACSN> <semi-colon separated list of PubMedID >
	columns are TAB separated
	each colum can contain SPACES

hugogenes.txt
	Positions of HUGO genes on Build37 (= hg19)
	From: http://chgr.mc.vanderbilt.edu/bushlab/wp-content/uploads/2011/06/hugogenes.txt
	<chromosome> <start position> <end position> <HUGO name>
	columns are TAB separated
	For chromosomes, 23 denotes chromosome X and 24 denotes chromosome Y

icogs_snp_list.csv
    From: http://ccge.medschl.cam.ac.uk/files/2014/03/icogs_snp_list.csv
