#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys
import numpy as np


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
