#ifndef maxflow_GT_h
#define maxflow_GT_h

#include "hi_pr.h"
#include "dimacs_parser.h"
#include "maxflow_solver.h"

class maxflow_GT : public maxflow_solver, public dimacs_parser_callback{
public:
	hi_pr g;
public:
       	maxflow_GT();
	void construct(const dynamic::num_array<int, 2> &E, dynamic::num_array<int, 2> &cap, dynamic::num_array<int, 1> &excess);
	tflow maxflow();
	void construct(const char *filename);
	virtual void save_cut_to_file(const std::string & filename)override;
	void save_cut();
	void save_nodes(int num_nodes = 0);
	virtual tflow cut_cost()override;
	virtual void get_cut(int * S)override;

public:
	virtual void allocate1(int n, int m, int S, int T) override;
	virtual void allocate2(int loop) override;
	virtual void read_arc(int loop, int u, int v, float cap1, float cap2) override;
};

#endif /* maxflow_GT_h */
