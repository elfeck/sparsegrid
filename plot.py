import numpy as np
import math
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from gridpoint import *

def gridpoint2d(point):
    fig = plt.figure()
    plt.plot(point.coord[0], point.coord[1], "ro")
    plt.axis([0, 1, 0, 1])
    plt.show()

gp = Gridpoint(2, [1,1], [1,1])
print(gp.coord)
gridpoint2d(gp)
