#!/usr/bin/env Rscript
# Visualize networks in dimacs format using igraph and save them in pdf

library("igraph")


#-------------------------------------------------
# get arg (file name and attraction) : 


#library("optparse")
# 
#option_list = list(
#    make_option(c("-f", "--fname"), type="character", default=NULL, 
#              help="dimacs fname name", metavar="character"),
#	 make_option(c("-a", "--attraction"), type="numeric", default=0.1, 
#              help="attraction of node [default= %default]", metavar="numeric")
#); 
# 
#opt_parser = OptionParser(option_list=option_list);
#args = parse_args(opt_parser);
#
#if (is.null(args$fname)){
#  print_help(opt_parser)
#  stop("At least one argument must be supplied (input dimacs fname).\n", call.=FALSE)
#}
# file = args$fname
# attraction=args$attraction


args = commandArgs(trailingOnly=TRUE)
# test if there is at least 1 argument: if not, return an error, else, take args
if (length(args)==0) {
    stop("dimacs fname is required.\n", call.=FALSE)
} else {    
    fname = args[1]    
    if (length(args)==1) {
        # default attration 
        args[2] = 0.1
    }
    attraction = args[2]    
}
out_fname = paste(fname,'.pdf', sep="") 

#---------------------
# Open the network 

g = read_graph(fname, format =  "dimacs", directed = F)

#---------------------
#Plot the network 

X11() # we want to see the plot even it is a script

#tkplot(g) 

l <- layout.drl(g, options=list(simmer.attraction=attraction))
plot(g, layout=l, vertex.size=3, vertex.label=NA)

dev.copy2pdf(file = out_fname) 
message (paste ("Network plot has been saved under", out_fname) )
#------------------
#wait
message("Press Return To Exit")
invisible(readLines("stdin", n=1))


