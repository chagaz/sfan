#ifndef dimacs_parser_h
#define dimacs_parser_h
#include "defs.h"

__forceinline int strtol10(char *& ps);

class dimacs_parser_callback
{
public:
	virtual void allocate1(int n ,int m, int S, int T)=0;
	virtual void read_arc(int loop, int u, int v, float cap1, float cap2)=0;
	virtual void allocate2(int loop){};
};

void dimacs_parser(const char * filename, dimacs_parser_callback & A, int loops = 1);

#endif
