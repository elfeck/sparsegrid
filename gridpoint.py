import numpy as np
import math
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class Gridpoint:

    def hat(l, i, x):
        return max(0, 1 - abs(2**l * x - i))

    def __init__(self, dim, level, index, funEval = None):
        self.dim = dim
        self.level = level
        self.index = index
        self.funEval = funEval

        self.alpha = 1
        self.tmpalpha = [0 for i in range(dim)]

        self.coord = []
        for i, l in zip(index, level):
            self.coord.append(float(i) / 2**l)


    def __str__(self):
        return self.toString()

    def toString(self):
        return str(self.level) + " | " + str(self.index)

    def getPhiAtTmp(self, dim, x):
        return (self.tmpalpha[dim-1] *
                Gridpoint.hat(self.level[dim-1], self.index[dim-1], x))

    def getPhi(self, dim, xval):
        return [Gridpoint.hat(self.level[dim - 1], self.index[dim - 1], x)
                for x in xval]

    def setFunEval(self, f):
        self.funEval = self.alpha = f(self.coord)
        self.tmpalpha = [self.alpha for i in self.tmpalpha]

    def adjTmpAlpha(self):
        sta = 0
        for ta in self.tmpalpha:
            sta += 0.5 * ta
        self.alpha = sta
        print("sta for " + str(self) + ": " + str(sta))

class GridHashMap:

    def hashGP(level, index):
        return "".join(list(map(str, level)) + list(map(str, index)))

    def unhashGP(string):
        return [list(map(int, string[0:2])), list(map(int, string[2:4]))]

    def __init__(self, maxlevel):
        self.maxlevel = maxlevel
        self.array = []
        self.hashMap = { }

    def __call__(self, level, index):
        return self.hashMap[GridHashMap.hashGP(level, index)]

    def insert(self, gp):
        self.array.append(gp)
        self.hashMap[GridHashMap.hashGP(gp.level, gp.index)] = gp

    def getSubspace(self, level):
        keys = list(self.hashMap.keys())
        return [self(GridHashMap.unhashGP(k)[0], GridHashMap.unhashGP(k)[1])
                for k in keys
                if k[0:2] == "".join(map(str, level))]

    def getChildren(self, gp, dim):
        level = list(gp.level)
        indices = list(gp.index)
        # new(i) = 2i - 1
        # new(i) = 2i + 1
        childIndices = [list(indices), list(indices)]
        childIndices[0][dim-1] = 2 * childIndices[0][dim-1] - 1
        childIndices[1][dim-1] = 2 * childIndices[1][dim-1] + 1
        level[dim-1] += 1
        keys = list(self.hashMap.keys())
        return [self(GridHashMap.unhashGP(k)[0], GridHashMap.unhashGP(k)[1])
                for k in keys
                if (k[0:2] == "".join(map(str, level)) and
                    (k[2:4] == "".join(map(str, childIndices[0])) or
                     k[2:4] == "".join(map(str, childIndices[1]))))]


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
            for x in point.getPhi(1, xval)
            for y in point.getPhi(2, yval)]
    if withAlpha:
        np.mult(point.alpha, zval)
    xmesh, ymesh = np.meshgrid(xval, yval)
    zmesh = np.reshape(zval, xmesh.shape)
    ax.plot_surface(xmesh, ymesh, zmesh, cstride=1, rstride=1)

def plot_subspaces2d(gridmap, withAlpha = False, step = 0.05):
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
                z1s = g.getPhi(1, xval)
                z2s = g.getPhi(2, yval)
                a = g.alpha if withAlpha else 1
                zval = np.add(zval, [a * z1 * z2 for z1 in z1s for z2 in z2s])
            zmesh = np.reshape(zval, xmesh.shape)
            ax.plot_surface(xmesh, ymesh, zmesh, cstride=1, rstride=1)
    return fig

def plot_summedSubspaces2d(gridmap, ax, withAlpha = False, step = 0.05):
    xval = yval = np.arange(0, 1 + step, step)
    xmesh, ymesh = np.meshgrid(xval, yval)
    zval = np.zeros(len(xval)**2)
    for g in gridmap.array:
        z1s = g.getPhi(1, xval)
        z2s = g.getPhi(2, yval)
        a = g.alpha if withAlpha else 1
        zval = np.add(zval, [a * z1 * z2 for z1 in z1s for z2 in z2s])
    zmesh = np.reshape(zval, xmesh.shape)
    ax.plot_surface(xmesh, ymesh, zmesh, cstride=1, rstride=1)

def plot_fun2d(f, ax, step = 0.05):
    xval = yval = np.arange(0, 1 + step, step)
    zval = [f([x,y])
            for x in xval
            for y in yval]
    xmesh, ymesh = np.meshgrid(xval, yval)
    zmesh = np.reshape(zval, xmesh.shape)
    ax.plot_surface(xmesh, ymesh, zmesh, cstride=1, rstride=1)

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

def buildSubspace2d(level, gridmap):
    pts_partial = [[i for i in range(1, 2**l) if i % 2 == 1] for l in level]
    pts = [[x, y] for x in pts_partial[0] for y in pts_partial[1]]
    for p in pts:
        gridmap.insert(Gridpoint(2, level, p))

def setAlphas(gridmap, f):
    # set function evaluations
    for gc in gridmap.array:
        gc.setFunEval(f)

    def setAlphaForGP(dim, par):
        print("Adjusting alpha for: " + str(par[-1]))
        children = gridmap.getChildren(par[-1], dim)
        for c in children:
            for p in par:
                adj = p.getPhiAtTmp(dim, c.coord[dim-1])
                c.tmpalpha[dim-1] -= adj
            setAlphaForGP(dim, list(par + [c]))

    def setForGP(dimIter, dimAdj, gp):
        print("Anchor for adj: " + str(gp))
        setAlphaForGP(dimAdj, [gp])
        children = gridmap.getChildren(gp, dimIter)
        for c in children:
            setForGP(dimIter, dimAdj, c)

    setForGP(1, 2, gridmap([1,1], [1,1]))
    setForGP(2, 1, gridmap([1,1], [1,1]))

    for gc in gridmap.array:
        gc.adjTmpAlpha()

def fun1(xv):
    x = xv[0]
    y = xv[1]
    mu = (1/2.0, 1/2.0)
    var = 0.02
    return math.exp(-(pow((x - mu[0]), 2) +
                      (pow((y - mu[1]), 2))) / (2 * var))

def fun2(xv):
    return (math.sin(2 * math.pi * xv[0]) * math.sin(math.pi * xv[1]))

fun = fun1

grid = fullgrid2d(3)
setAlphas(grid, fun)
fig1 = plot_subspaces2d(grid, True)
fig1.show()

fig2 = plt.figure()
ax21 = fig2.add_subplot(111, projection="3d")
plot_fun2d(fun, ax21)
fig2.show()

fig3 = plt.figure()
ax31 = fig3.add_subplot(111, projection="3d")
plot_summedSubspaces2d(grid, ax31, True)
fig3.show()
