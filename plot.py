#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys
import numpy as np
from glob import glob

def plot_edges(filename, figure=None, point_style='k.', edge_style='k-'):
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
    if figure is None:
        # Create the figure
        p = plt.figure()
    else:
        p = figure
    plt.plot(verts[:,0], verts[:,1], point_style, markersize=15)
    plt.plot(edgex.T, edgey.T, edge_style)
    return p

def plot_edges_parallel(prefix, point_style=('k.',), edge_style=('k-',)):
    p = plt.figure()
    for i, f in enumerate(glob(prefix + "_*.txt")):
        ps = point_style[i] if len(point_style) > 1 else point_style[0]
        es = edge_style[i] if len(edge_style) > 1 else edge_style[0]
        plot_edges(f, figure=p, point_style=ps, edge_style=es)
    return p

def read_points(file_name, delim=' '):
    with open(file_name, mode='r') as f:
        # Skip the header
        f.readline()
        return np.asarray([(float(l.split(delim)[0]), float(l.split(delim)[1])) for l in f])

def plot_partition1d(prefix, partition, delim=' ', ext='txt', point_style='k.', line_style='r-', **kwargs):
    point = np.concatenate([read_points(f, delim=delim) for f in glob(prefix+"_*."+ext)], axis=0)
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


def plot_points(filename, delim=' ', point_style='b.', **kwargs):
    p = plt.figure()
    point = read_points(filename, delim=delim)
    plt.plot(point[:,0], point[:,1], point_style)
    if 'xlim' in kwargs:
        plt.xlim(kwargs['xlim'])
    if 'ylim' in kwargs:
        plt.xlim(kwargs['ylim'])
    return p

def plot_partitioned_points(prefix, delim=' ', ext='txt', point_styles=('bs', 'r^', 'gv', 'cx'), **kwargs):
    p = plt.figure()
    for i, f in enumerate(glob(prefix+"_*."+ext)):
        plot_points(f, delim=delim, point_style=point_styles[i], **kwargs)
    return p
