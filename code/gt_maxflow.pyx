cdef extern from "c++sources/entry_point.cpp":
	void entry_point(int argc, char *argv, int num_nodes)
    
cpdef python_entry_point(char *argv, int num_nodes):
	entry_point(1, argv, num_nodes)