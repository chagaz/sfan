#include "gt_maxflow_sources/GT.cpp"

int entry_point(int argc, char *argv, int num_nodes);

int entry_point(int argc, char *argv, int num_nodes)
{	
    int arguments_count = argc + 1;
    const char *arguments_value[] = {"gt_maxflow\0", static_cast<const char*>(argv)};
    
    return (gt5_main(arguments_count, arguments_value, num_nodes));
}
