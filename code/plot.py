import matplotlib.pyplot as plt
import numpy as np
import random

# known bug : 
# should not be used with pdb debugger...

# TODO : fill box + legende
# TODO : add means as star
# TODO : add all points 
# TODO : integrer plot.py code 


algos_names = ('SConES', 'MSConESnp', 'MSConES')

colors = ['darkkhaki', 'royalblue', 'white']




def gene_fake_data (num_tasks) : 
    # Generate data
    data = {}
    for key in algos_names : 
        data[key] = []
    n = 5
    for task_id in xrange(num_tasks) :
        upper = random.randint(0, 1000)  
        data[algos_names[0]].append(np.random.uniform(0, upper, size=n))
        data[algos_names[1]].append(np.random.uniform(0, upper, size=n-1))
        data[algos_names[2]].append(np.random.uniform(0, upper, size=n+103))
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
    


def vertical_boxplots(data, name) : 
    num_tasks = len(data['MSConES'] ) 

    fig, axes = plt.subplots(ncols=num_tasks, sharey=True)
    fig.subplots_adjust(wspace=0)
    fig.canvas.set_window_title(name)

    for i, (ax, task_id) in enumerate( zip(axes, xrange(num_tasks) )) : 
        print task_id
        # plot task per task : 
        boxp = ax.boxplot([data[algo_id][task_id] for algo_id in algos_names], vert = True)
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

    axes[num_tasks / 2 -1].set_title(name)
    #axes[num_tasks / 2].set_xlabel('Task')
    axes[0].set_ylabel('Measure')

    #make_lined_legend()
    #make_filled_legend() 
    
    fig.tight_layout() #ajuste le cadrage
    plt.savefig(name+'.png')   
    #plt.show()



def horizontal_boxplots(data,name) : 
    num_tasks = len(data['MSConES'] ) 

    fig, axes = plt.subplots(num_tasks, 1, sharex=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.canvas.set_window_title(name)
    
    for i, (ax, task_id) in enumerate( zip(axes, xrange(num_tasks) )) : 
        print task_id
        # plot task per task : 
        boxp = ax.boxplot([data[algo_id][task_id] for algo_id in algos_names], vert = False)
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

    axes[0].set_title(name)
    #axes[num_tasks / 2].set_xlabel('Task')
    axes[num_tasks / 2 -1].set_ylabel('Measure')

    #make_lined_legend()
    #make_filled_legend() 
    
    fig.tight_layout() #ajuste le cadrage
    plt.savefig(name+'.png')   
    #plt.show()


def vertical_barplots(data, name) : 
    num_tasks = len(data['MSConES'] ) 

    fig, axes = plt.subplots(ncols=num_tasks, sharey=True)
    fig.subplots_adjust(wspace=0)
    fig.canvas.set_window_title(name)

    x_loc = np.arange( len(algos_names) ) # [0, 1, 2] # 3 algos
    bar_width = 0.35

    for i, (ax, task_id) in enumerate( zip(axes, xrange(num_tasks) )) : 
        ax.bar(x_loc, [np.mean(data[algo_id][task_id]) for algo_id in algos_names], bar_width )
        ax.set_xticks(x_loc)
        ax.set(xticklabels=algos_names, xlabel=task_id)
        #ax.margins(0.05) # Optional

    fig.tight_layout() #ajuste le cadrage
    plt.savefig(name+'.png')    
    #plt.show()





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

################################################################################
################################################################################


if __name__ == "__main__":

    num_tasks = 4


    # data [algo_id][task_id] = list of values
    data = gene_fake_data(num_tasks)

    # create fig instances : 

    #-----------
    # one for boxplots : 
    # ... VERTICAL boxplots : 
    vertical_boxplots(data, 'vertical_boxplots') 
    # ... HORIZONTAL boxplots  : 
    horizontal_boxplots(data, 'horizontal_boxplots') 
    #-----------
    # one for barplots : 
    vertical_barplots(data, 'vertical_barplots') 
    
    print("plot, THE END")

