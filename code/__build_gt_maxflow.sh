#!/bin/bash

rm -rf ./gt_maxflow.so

cython --cplus gt_maxflow.pyx # generate gt_maxflow.cpp

g++ -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 \
	-o gt_maxflow.so gt_maxflow.cpp \
\
	c++sources/gt_maxflow_sources/block_allocator.cpp \
	c++sources/gt_maxflow_sources/num_array.cpp \
	c++sources/gt_maxflow_sources/logs.cpp \
	c++sources/gt_maxflow_sources/performance.cpp \
	c++sources/gt_maxflow_sources/pvect.cpp \
	c++sources/gt_maxflow_sources/xfs.cpp \
	c++sources/gt_maxflow_sources/xstringstream.cpp \
	c++sources/gt_maxflow_sources/file_stream.cpp \
	c++sources/gt_maxflow_sources/binary_stream.cpp \
	c++sources/gt_maxflow_sources/text_stream.cpp \
\
	c++sources/gt_maxflow_sources/dimacs_parser.cpp \
	c++sources/gt_maxflow_sources/maxflow_solver.cpp \
\
	c++sources/gt_maxflow_sources/maxflow_GT.cpp \
	c++sources/gt_maxflow_sources/construct.cpp \
	c++sources/gt_maxflow_sources/hi_pr.cpp \
	c++sources/gt_maxflow_sources/parser.cpp \
	c++sources/gt_maxflow_sources/timer.cpp \
\
	c++sources/gt_maxflow_sources/SFILE.cpp \
\
	-w -O2 -fpermissive -fPIC -fomit-frame-pointer -Winline --param inline-unit-growth=1000 --param max-inline-insns-single=1000 --param large-function-growth=500 -DMX_COMPAT_32

rm -rf ./build ./gt_maxflow.cpp
