#include "maxflow_solver.h"
#include "maxflow_GT.h"
#include "xfs.h"

using namespace debug;
using namespace dynamic;

debug::PerformanceCounter c1;

void solve(const char * file){
    try{
	maxflow_GT * solver;
	solver = 0;
	dimacs_parser_callback * constructor;
	solver = new maxflow_GT;
	solver -> g.globUpdtFreq = 0.5;
	constructor = solver;
	
	// Parse .dimacs
	PerformanceCounter c0;
	c0.start();
	dimacs_parser(file,*constructor,2);
	c0.stop();
	debug::stream << "parser: " << c0.time() <<"\n";
	
	// Compute max flow
	double F = solver->maxflow();
	debug::stream << "F: " << F << "\n"; 
	solver->print_info();
	txt::StringStream log_name;
	log_name << std::string(file) << ".GT05.cut";
	solver->save_cut_to_file(log_name);
	
	//cleanup
	delete solver;
    } catch(const std::exception & e){
	debug::stream<<"Exception: "<< e.what() <<"\n";
	exit(1);
    } catch(...){
	debug::stream<<"ERROR\n";
	exit(1);
    };
}

void try_solve(const char * file){
    solve(file);
};

int main(int argc, char *argv[]){
    const char * solver;
    const char * file;
    const char * options;
    std::string root;
    dynamic::fixed_array1<std::string> solvers;
    solvers.reserve(20);
    root = xfs::getPath(argv[0]);
    printf("my path: %s\n",root.c_str());
    file =  "test/BVZ-tsukuba0.max";
    solvers.push_back("GT05");
    
    txt::StringStream log_name;
    log_name << std::string(file) << ".GT05.sol";
    txt::FileStream f(log_name.c_str());
    txt::EchoStream g(&f,&stdoutStream::the_stream);
    debug::stream.attach(&g);
    try_solve(file);
    debug::stream<<"All Ok\n";
    
    return 0;
};
