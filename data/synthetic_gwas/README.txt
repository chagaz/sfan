FILE FORMAT

.gmt:
	Each line corresponds to one gene
	<gene name in ACSN> <dataset name> <gene name 1 in HUGO>  <gene name 2 in HUGO> ...
	columns are TAB separated
	each colum can contain SPACES

.sif:
	Each line corresponds to one interaction (network edge)
	<gene 1 name in ACSN> <name of relationship> <gene 2 name in ACSN>
	columns are TAB separated
	each colum can contain SPACES


hugogenes.txt
	Positions of HUGO genes on Build37 (= hg19)
	From: http://chgr.mc.vanderbilt.edu/bushlab/wp-content/uploads/2011/06/hugogenes.txt
	<chromosome> <start position> <end position> <HUGO name>
	columns are TAB separated
	For chromosomes, 23 denotes chromosome X and 24 denotes chromosome Y
