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

mart_export.txt
	Human genes list and positions with HUGO genes symbol on Build37 (= hg19)
	From: http://grch37.ensembl.org/biomart/martview/
	Dowloaded: 11 june 2016
	Ensembl GRCh37 Biomart Tool 
        Database : Ensembl Genes 
        Dataset : Homo sapiens genes (GRCh37.p13) 
        Attributes : 
	        Gene > Chromosome Name, Gene Start (bp), Gene End (bp) # The entire gene and not only the exons #
	        External : HGNC Symbol
        Unique Results only : yes
	<chromosome> <start position> <end position> <HUGO symbol>
	columns are TAB separated

icogs_snp_list.csv
    From: http://ccge.medschl.cam.ac.uk/files/2014/03/icogs_snp_list.csv
