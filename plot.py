#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys
import numpy as np
from glob import glob

def plot_edges(filename, point_style='k.', edge_style='k-'):
    # Read file
    with open(filename, mode='r') as f:
        ndim = 2
        nVerts, nEdges = (int(x) for x in f.readline().split())
        verts = np.empty((nVerts, ndim))
        for i in range(nVerts):
            verts[i,:] = [float(x) for x in f.readline().split()]
            edges = [0] * nEdges
        for i in range(nEdges):
            edges[i] = [int(x) for x in f.readline().split()]
    # Edge arrays for plotting
    edgex = np.asarray([[verts[e[0],0], verts[e[1],0]] for e in edges])
    edgey = np.asarray([[verts[e[0],1], verts[e[1],1]] for e in edges])
    # Create the figure
    p = plt.figure()
    plt.plot(verts[:,0], verts[:,1], point_style, markersize=15)
    plt.plot(edgex.T, edgey.T, edge_style)
    return p

def _read_csv(file_name):
    with open(file_name, mode='r') as f:
        # Skip the header
        f.readline()
        return np.asarray([(float(l.split(',')[0]), float(l.split(',')[1])) for l in f])

def plot_partition1d(prefix, partition, point_style='k.', line_style='r-', **kwargs):
    point = np.concatenate([_read_csv(f) for f in glob(prefix+"_*.csv")], axis=0)
    xmin = point[:,0].min()
    xmax = point[:,0].max()
    ymin = point[:,1].min()
    ymax = point[:,1].max()
    p = plt.figure()
    plt.plot(point[:,0], point[:,1], point_style)
    if 'xlim' in kwargs:
        plt.xlim(kwargs['xlim'])
        xmin = kwargs['xlim'][0]
        xmax = kwargs['xlim'][1]
    if 'ylim' in kwargs:
        plt.xlim(kwargs['ylim'])
        ymin = kwargs['ylim'][0]
        ymax = kwargs['ylim'][1]

    for bound in partition:
        plt.plot([bound, bound], [ymin, ymax], line_style)
    return p

def plot_partitioned_points(prefix, point_styles=('bs', 'r^', 'gv', 'cx'), **kwargs):
    p = plt.figure()
    for i, f in enumerate(glob(prefix+"_*.csv")):
        point = _read_csv(f)
        plt.plot(point[:,0], point[:,1], point_styles[i])
    if 'xlim' in kwargs:
        plt.xlim(kwargs['xlim'])
    if 'ylim' in kwargs:
        plt.xlim(kwargs['ylim'])
    return p
