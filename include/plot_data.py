#!/usr/bin/env python

from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


def plot_x_ny(x, y, x_label, y_label):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    for y_axis in y:
        ax.plot(x, y_axis, '-')

    #ax.set_yscale('log')
    #plt.tight_layout()
    plt.show()


