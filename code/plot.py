import matplotlib.pyplot as plt
import numpy as np

def bar_plot(measure, f_name ):
    # code inspired of : http://matplotlib.org/examples/api/barchart_demo.html
    

    plt.style.use('ggplot')
    
    plt.figure()
    algos = ['st', 'np', 'msfan']
    means_to_plot = {}
    std_to_plot = {}
    with open(f_name, 'r') as f :
        for i, line in enumerate(f) :
            content = line.split('|') 
            means_to_plot[algos[i]] = [float(item) for item in content[0].split() ]
            std_to_plot[algos[i]] = [float (item) for item in content[1].split() ]
    
    num_groups = len(means_to_plot[algos[i]]) # = num_tasks
    labels_groups = [str(group) for group in xrange(num_groups) ]
    index = np.arange(num_groups)  # the x locations for the groups

    bar_width = 0.2    # the width of the bars <= 1 / float ( len(algos) )
    colors = ['b', 'r', 'k'] # blue, red, black
    opacity = 0.5 #0.4

    error_config = {'ecolor': '0.3'} #TODO : kezako ? 


    plt.xlabel('Task')
    plt.ylabel(measure)
    plt.xticks(index + bar_width/2 + bar_width, labels_groups)
    
    plt.title('%s by task and algo' % measure)

    for i, algo in enumerate(algos) : 
        # compute bars for each algos 
        print i, algo 
        print means_to_plot[algo]
        print std_to_plot[algo]
        plt.bar (index + (i%len(algos))*bar_width, 
                means_to_plot[algo], 
                bar_width,
                alpha=opacity,
                color=colors[i%len(algos)],
                yerr=std_to_plot[algo],
                error_kw=error_config,
                label=algo
        )

    lgd = plt.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))
    plt.tight_layout()
    
    plt.savefig('plot%s.png'%measure, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')


if __name__ == "__main__":
    #TODO : check arg
    bar_plot()
    print("THE END")
        


