#include "dimacs_parser.h"
#include <stdio.h>
#include "dynamic.h"
#include <string.h>

#include "SFILE.hpp"
using namespace KCS;

class buff_read{
public:
    static const int			 max_size = 1024*1024*32;	//8Mb
    dynamic::fixed_array1<char>	 buff;
    int				 sz;
    int				 pos;
    SFILE				*f;
private:
    void readahead(){
	sz					  = f->read(buff.begin(), 1, buff.size());
	pos					  = 0;
    };
    __forceinline bool ready(){
	if(pos<sz)return true;
	
	if (f->eos())
	    return false;
	
	readahead();
	return true;
    };
public:
    bool eof() const
    {
	return (pos == sz && f->eos());
    }
    buff_read(SFILE *_f) : 
	f(_f)
    {
	pos = 0;
	sz = 0;
	buff.resize(max_size);
    };
    //
    ~buff_read(){
    };
    __forceinline char getc(){
	ready();
	return buff[pos++];
    };
    //
    __forceinline void skip_line(){
	while(ready()){
	    while(pos<sz){
		if(buff[pos]=='\n')return;
		++pos;
	    };
	};
    };
    //
    __forceinline int read_line(char * s,int max_count){
	int k = 0;
	while(ready() && k<max_count-1){
	    while(pos<sz && k<max_count-1){
		s[k] = buff[pos];
		if(buff[pos]=='\n'){
		    s[k+1] = 0;
		    return k;
		};
		++pos;
		++k;
	    };
	};
	s[k+1] = 0;
	return k;
    };
};

__forceinline char * eatspace(char * ps){
    do{
	char c = *ps;
	if(c==' ' || c== '\t' || c== '\n' || c=='\r'){
	    ++ps;
	    continue;
	};
	return ps;
    }while(1);
};
__forceinline int strtol10(char *& ps){
    int r = 0;
    bool isnumber = false;
    ps = eatspace(ps);
    do{
	char c = *ps;
	if(c<'0' || c>'9'){
	    if(!isnumber)
		std::cout << "number expected" << std::endl;
	    return r;
	};
	isnumber = true;
	r = r*10+(c-'0');
	++ps;
    }while(1);
};

void dimacs_parser(const char * filename, dimacs_parser_callback & A, int loops)
{
    SFILE f(filename);
    
    if (f.operator!())
	std::cout << "cant read c_string" << std::endl;
    
    char s[1024];
    char s1[1024];
    long int pos;// remember position where arcs start in the stream
    char c = f.getc();
    int n, m;
    int S, T;
    struct tlast_arc
    {
	int u;
	int v;
	float cap;
	bool parsed;
	tlast_arc():parsed(true){};
    } last_arc;
    
    for ( ; !f.eos() && c != 'a'; c = f.getc())
	{
	    switch (c)
		{
		case 'p'://problem specification
		    f.scans(" %*s %i %i ", &n, &m); // number of nodes and arcs
		    break ;
		case 'c'://comment line, look for regulargreed, complexgrid
		    f.gets(s, 1024);//read comment line until the end
		    sscanf(s, " %s ", s1);
		    break ;
		case 'n'://source or sink nodes
		    int v;
		    f.scans(" %i %c ", &v, s);
		    --v; //zero-based index
		    if (s[0] == 's')
			S = v;
		    else
			T = v;
		    break ;
		}
	}
    pos = f.tell() - 1;
    A.allocate1(n, m, S, T);

    // now a double loop over edges - count and read
    for (int loop = 0; loop < loops; ++loop)
	{
	    f.seek(pos, SEEK_SET); //rewind to where arcs begin
	    buff_read ff(&f);
	    while (!ff.eof())
		{
		    c = ff.getc();
		    switch (c)
			{
			case 'c':
			    ff.skip_line();
			    break ;
			case 'a':
			    int head;
			    int tail;
			    float cap;
			    ff.read_line(s, 1024);
			    char *ps = s;
			    head = strtol10(ps);
			    tail = strtol10(ps);
			    cap = atof(ps);
			    if (cap == 0.)
				{
				    continue ;//skip zero arcs
				}
			    --head;//to zero-based index
			    --tail;
			    //break;
			    //if u is source or v is sink, assign this capacity to the excess
			    //A.read_arc(loop,u,v,cap,0);//next arc was a reverse one
			    //continue;
			    if (head == S || tail == T)
				A.read_arc(loop, head, tail, cap, 0);
			    else
				{
				    if (!last_arc.parsed)
					{
					    if (last_arc.u == tail && last_arc.v == head)
						{
						    A.read_arc(loop, head, tail, cap, last_arc.cap);//next arc was a reverse one
						    last_arc.parsed = true;
						}
					    else
						{
						    A.read_arc(loop, last_arc.u, last_arc.v, last_arc.cap, 0);//parse as unpaired
						    last_arc.u = head;
						    last_arc.v = tail;
						    last_arc.cap = cap;
						    last_arc.parsed = false;
						}
					}
				    else
					{
					    last_arc.u = head;
					    last_arc.v = tail;
					    last_arc.cap = cap;
					    last_arc.parsed = false;
					}
				}
			    break ;//case
			}
		}
	    if(!last_arc.parsed)
		{
		    A.read_arc(loop, last_arc.u, last_arc.v, last_arc.cap, 0);//parse as unpaired
		    last_arc.parsed = true;
		}
	    A.allocate2(loop);
	}
    f.close();
}
