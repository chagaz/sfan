import matplotlib.pyplot as plt
import numpy as np
import random


algos_names = ('SConES', 'MSConESnp', 'MSConES')
num_tasks = 4
colors = ['darkkhaki', 'royalblue', 'white']



def vertical_boxplots(data) : 

    fig, axes = plt.subplots(ncols=num_tasks, sharey=True)
    fig.subplots_adjust(wspace=0)
    fig.canvas.set_window_title('Boxlots')

    for i, (ax, task_id) in enumerate( zip(axes, xrange(num_tasks) )) : 
        print task_id
        # plot task per task : 
        boxp = ax.boxplot([data[task_id][algo_id] for algo_id in algos_names], vert = True)
        # add vertical x ticks : 
        ax.set(xticklabels=algos_names)#, xlabel=task_id)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
        ax.margins(0.05) # Optional
        # set colors : 
        plt.setp(boxp['boxes'], color='black')
        plt.setp(boxp['whiskers'], color='black')
        plt.setp(boxp['fliers'], color='red', marker='+')
        # set grid : 
        add_hlines(ax) 

    axes[num_tasks / 2 -1].set_title('Title')
    #axes[num_tasks / 2].set_xlabel('Task')
    axes[0].set_ylabel('Measure')

    #make_lined_legend()
    #make_filled_legend() 
    
    fig.tight_layout() #ajuste le cadrage
    plt.savefig('boxplots1.png')   
    plt.show()



def horizontal_boxplots(data) : 
    fig, axes = plt.subplots(num_tasks, 1, sharex=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.canvas.set_window_title('Boxlots')

    for i, (ax, task_id) in enumerate( zip(axes, xrange(num_tasks) )) : 
        print task_id
        # plot task per task : 
        boxp = ax.boxplot([data[task_id][algo_id] for algo_id in algos_names], vert = False)
        # left x labels  = algos names  
        ax.set(yticklabels=algos_names)
        # right x labels = num task : 
        rax =  ax.twinx()
        rax.set(yticklabels=[])
        rax.set_ylabel(task_id, rotation=0)
        #
        ax.margins(0.05) # Optional
        # set colors : 
        plt.setp(boxp['boxes'], color='black')
        plt.setp(boxp['whiskers'], color='black')
        plt.setp(boxp['fliers'], color='red', marker='+')
        # set grid : 
        add_vlines(ax) 

    axes[0].set_title('Title')
    #axes[num_tasks / 2].set_xlabel('Task')
    axes[3].set_xlabel('Measure')

    #make_lined_legend()
    #make_filled_legend() 
    
    fig.tight_layout() #ajuste le cadrage
    plt.savefig('boxplots2.png')   
    plt.show()

