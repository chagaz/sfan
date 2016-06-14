import matplotlib.pyplot as plt
import numpy as np
import random


algos_names = ('SConES', 'MSConESnp', 'MSConES')
num_tasks = 4
colors = ['darkkhaki', 'royalblue', 'white']




def gene_fake_data () : 
    # Generate data
    data = {}
    for key in xrange(num_tasks) : 
        data[key] = {}
    n = 5
    for k,v in data.iteritems():
        upper = random.randint(0, 1000)
        v[algos_names[0]] = np.random.uniform(0, upper, size=n)
        v[algos_names[1]] = np.random.uniform(0, upper, size=n-1)
        v[algos_names[2]] = np.random.uniform(0, upper, size=n+103)
    return data
####################################################




def add_hlines(axe) : 
    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    axe.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)
    # Hide these grid behind plot objects
    axe.set_axisbelow(True)

def add_vlines(axe) : 
    # Add a vertical grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    axe.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)
    # Hide these grid behind plot objects
    axe.set_axisbelow(True)
    


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


def vertical_barplots(data) : 
    fig, axes = plt.subplots(ncols=num_tasks, sharey=True)
    fig.subplots_adjust(wspace=0)
    fig.canvas.set_window_title('Barplots')

    x_loc = np.arange( len(algos_names) ) # [0, 1, 2] # 3 algos
    bar_width = 0.35

    for i, (ax, task_id) in enumerate( zip(axes, xrange(num_tasks) )) : 
        ax.bar(x_loc, [np.mean(data[task_id][algo_id]) for algo_id in algos_names], bar_width )
        ax.set_xticks(x_loc)
        ax.set(xticklabels=algos_names, xlabel=task_id)
        #ax.margins(0.05) # Optional

    fig.tight_layout() #ajuste le cadrage
    plt.savefig('barplots.png')    
    plt.show()





def make_lined_legend():
    # draw temporary red and blue lines and use them to create a legend
    hM, = plt.plot([1,1],'m-') #magenta
    hB, = plt.plot([1,1],'b-') #blue
    hK, = plt.plot([1,1],'k-') #black
    plt.legend((hM, hB, hK),algos_names)
    # on ne veut pas vraiment les dessiner : 
    hO.set_visible(False) 
    hP.set_visible(False)
    hB.set_visible(False)


def make_filled_legend():
    plt.figtext(0.80, 0.18,  algos_names[0],
                backgroundcolor=colors[0], color='black', weight='roman',
                size='x-small')
    plt.figtext(0.80, 0.145, algos_names[1],
                backgroundcolor=colors[1], color='white', weight='roman', 
                size='x-small')
    plt.figtext(0.80, 0.110, algos_names[2],
                 backgroundcolor=colors[2], color='black', weight='roman', 
                size='x-small')
    plt.figtext(0.80, 0.015, '*', color='white', backgroundcolor='silver',
                weight='roman', size='medium')
    plt.figtext(0.83, 0.015, 'Mean', color='black', weight='roman',
                size='x-small')
            
