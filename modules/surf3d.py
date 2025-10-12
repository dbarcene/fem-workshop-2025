#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : surf3d.py
# Author            : David Barcene <david.barcene@utp.ac.pa>
# Date              : 18.07.2022
# Last Modified Date: 18.07.2022
# Last Modified By  : David Barcene <david.barcene@utp.ac.pa>


from matplotlib.figure import projections
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def surf3d(xx, yy, V):

    fig = plt.figure(figsize=(14, 9))
    ax = plt.axes(projection='3d')
    ax.plot_surface(xx, yy, V, cmap=plt.cm.plasma)
    ax.set_zlim(0, 5e-5)
    plt.show()
