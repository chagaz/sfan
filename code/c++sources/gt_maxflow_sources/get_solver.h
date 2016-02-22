#include "block_allocator.h" //memory manager initializes first
#include "region_splitter.h"
#include "maxflow_solver.h"

void get_solver(const char * _solver, region_graph * G, region_splitter * splitter, maxflow_solver *& solver, dimacs_parser_callback *& constructor, intf & slice, bool unload, const char * tmp="");