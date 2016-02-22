#include "maxflow_solver.h"
#include "maxflow_GT.h"
#include "xfs.h"

#include <iostream>
#include <fstream>

using namespace debug;

debug::PerformanceCounter c1;

#include <sys/time.h>


void solve(const char * file, int num_nodes){
    try{
	
	std::streambuf* cout_sbuf = std::cout.rdbuf();
	std::ofstream   fout("/dev/null");
	std::cout.rdbuf(fout.rdbuf());
	
	maxflow_GT * solver;
	dimacs_parser_callback * constructor;

	struct timeval tbeginalgo, tendalgo, tbeginread, tendread;
	double texecalgo(0), texecread(0);

	solver = new maxflow_GT;
	((maxflow_GT*)solver)->g.globUpdtFreq = 0.5;
	constructor = (maxflow_GT*)solver;

	gettimeofday(&tbeginread, NULL);
	dimacs_parser(file,*constructor,2);
	gettimeofday(&tendread, NULL);
	texecread = (double)(1000 * (tendread.tv_sec - tbeginread.tv_sec) + ((tendread.tv_usec - tbeginread.tv_usec)/1000));
	debug::stream << "exec time for building graph = " << texecread << " ms\n";

	gettimeofday(&tbeginalgo, NULL);    
	double F = solver->maxflow();
	gettimeofday(&tendalgo, NULL);
	texecalgo = (double)(1000 * (tendalgo.tv_sec - tbeginalgo.tv_sec) + ((tendalgo.tv_usec - tbeginalgo.tv_usec)/1000));
	debug::stream << "exec time for maxflow = " << texecalgo << " ms\n";
	solver->print_info(); 

	solver->save_cut();
	
	std::cout.rdbuf(cout_sbuf);
	solver->save_nodes(num_nodes);
    
	delete solver;
    } catch(...){
	debug::stream << "ERROR\n";
	exit(1);
    };
}

int gt5_main(int argc, char *argv[], int num_nodes){
    if (argc != 2)
	{
	    std::cout << "ERROR Wrong number of arguments in main" << std::endl;
	}
    else
	{
	    const char * file;
	    file = argv[1];
	    solve(file, num_nodes);
      
	    return 0;
	};
}
