import numpy as np
import math
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class Gridpoint:

    def hat(l, i, x):
        return max(0, 1 - abs(2**l * x- i))

    def __init__(self, dim, levels, indices, funEval = None):
        self.dim = dim
        self.levels = levels
        self.indices = indices
        self.funEval = funEval

        self.alpha = 1

        self.coord = []
        for i, l in zip(indices, levels):
            self.coord.append(float(i) / 2**l)

    def evalAt(self, x):
        result = self.alpha
        for d in range(self.dim):
            result *= Gridpoint.hat(self.levels[i], self.indices[i], x[i])
        return result

    def getPhi(self, dim, xval, withAlpha = False):
        a = self.alpha if withAlpha else 1.0
        return [a*Gridpoint.hat(self.levels[dim - 1], self.indices[dim - 1], x)
                for x in xval]

    def setFunEval(self, f):
        self.funEval = f(self.coord)

class GridHashMap:

    def hashGP(levels, indices):
        return "".join(list(map(str, levels)) + list(map(str, indices)))

    def unhashGP(string):
        return [list(map(int, string[0:2])), list(map(int, string[2:4]))]

    def __init__(self, maxlevel):
        self.maxlevel = maxlevel
        self.array = []
        self.hashMap = { }

    def __call__(self, levels, indices):
        return self.hashMap[GridHashMap.hashGP(levels, indices)]

    def insert(self, gp):
        self.array.append(gp)
        self.hashMap[GridHashMap.hashGP(gp.levels, gp.indices)] = gp

    def getSubspace(self, level):
        keys = list(self.hashMap.keys())
        return [self(GridHashMap.unhashGP(k)[0], GridHashMap.unhashGP(k)[1])
                for k in keys
                if k[0:2] == "".join(map(str, level))]


def plot_gridpoint2d(point, ax):
    ax.plot(point.coord[0], point.coord[1], "ro")
    ax.axis([0, 1, 0, 1])

def plot_phi1d(point, dim, ax, step = 0.01):
    xval = np.arange(0, 1 + step, step)
    yval = point.getPhi(dim, xval)
    ax.plot(xval, yval)
    ax.axis([0, 1, 0, 1])

def plot_phi2d(point, ax, withAlpha = False, step = 0.05):
    xval = yval = np.arange(0, 1 + step, step)
    zval = [x * y
            for x in point.getPhi(1, xval, withAlpha)
            for y in point.getPhi(2, yval, withAlpha)]
    xmesh, ymesh = np.meshgrid(xval, yval)
    zmesh = np.reshape(zval, xmesh.shape)
    ax.plot_surface(xmesh, ymesh, zmesh, cstride=1, rstride=1)

def plot_subspaces(gridmap, withAlpha = False, step = 0.05):
    fig = plt.figure()
    mlx = gridmap.maxlevel[0]
    mly = gridmap.maxlevel[1]

    xval = yval = np.arange(0, 1 + step, step)
    xmesh, ymesh = np.meshgrid(xval, yval)
    for lx in range(1, mlx + 1):
        for ly in range(1, mly + 1):
            ax = fig.add_subplot(mlx, mly, (ly - 1) * mlx + lx,
                                 projection="3d")
            subspace = gridmap.getSubspace([lx, ly])
            zval = np.zeros(len(xval)**2)
            for g in subspace:
                z1s = g.getPhi(1, xval, withAlpha)
                z2s = g.getPhi(2, yval, withAlpha)
                zval = np.add(zval, [z1 * z2 for z1 in z1s for z2 in z2s])
            zmesh = np.reshape(zval, xmesh.shape)
            ax.plot_surface(xmesh, ymesh, zmesh, cstride=1, rstride=1)
    return fig

def fullgrid2d(maxlevel):
    gridmap = GridHashMap([maxlevel, maxlevel])
    for lx in range(1, maxlevel + 1):
        for ly in range(1, maxlevel + 1):
            buildSubspace2d([lx, ly], gridmap)
    return gridmap

def sparsegrid2d(maxlevel):
    gridmap = GridHashMap([maxlevel, maxlevel])
    for lx in range(1, maxlevel + 1):
        for ly in range(1, maxlevel + 1):
            if lx + ly > maxlevel + 1:
                continue
            buildSubspace2d([lx, ly], gridmap)
    return gridmap

def buildSubspace2d(levels, gridmap):
    pts_partial = [[i for i in range(1, 2**l) if i % 2 == 1] for l in levels]
    pts = [[x, y] for x in pts_partial[0] for y in pts_partial[1]]
    for p in pts:
        gridmap.insert(Gridpoint(2, levels, p))

def fun(xv):
    x = xv[0]
    y = xv[1]
    mu = (1/2.0, 1/2.0)
    var = 0.02
    return math.exp(-(pow((x - mu[0]), 2) +
                      (pow((y - mu[1]), 2))) / (2 * var))


grid = sparsegrid2d(3)
fig = plot_subspaces(grid)
fig.show()
