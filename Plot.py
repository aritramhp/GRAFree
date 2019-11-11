#!/usr/bin/env python
# ========
# PLOT
# ========

import Header
from Header import *


# -------------------------------------------------------------------------------------
# This function call different plotting functions for different type of plotting
# -------------------------------------------------------------------------------------
def mainPlotting(x_axis, y_axis, graph_filename, block_size=0):
    # collect x axis data
    value = x_axis
    val = map(float, value)
    x_axis = val
    # collect y axis data
    value = y_axis
    val = map(float, value)
    y_axis = val
    
    # Plot in graph
    graph_dirname = graph_filename.split('/')[:-1]
    graph_dirname = '/'.join(graph_dirname)
    
    if not os.path.isdir(graph_dirname):
        mkdr_cmd = 'mkdir -p ' + graph_dirname
        os.system(mkdr_cmd)
    
    x_axis = numpy.array(x_axis)
    y_axis = numpy.array(y_axis)
    
    plt.rcParams.update({'font.size': 18})
    plt.gca().set_color_cycle(['blue'])
    plt.plot(x_axis, y_axis, linewidth=.25)
    plt.plot(x_axis[0], y_axis[0], 'og')
    plt.plot(x_axis[-1], y_axis[-1], 'or')
    plt.grid()
    # plt.show()
    if not block_size == 0:
        plt.xlim(-1 * block_size, block_size)
        plt.ylim(-1 * block_size, block_size)
    plt.savefig(graph_filename)
    plt.close()
    
    del x_axis, y_axis
    gc.collect()
