# libcgeom

[![Build Status](https://travis-ci.org/jdahm/libcgeom.svg?branch=master)](https://travis-ci.org/jdahm/libcgeom)

These aren't the droids you seek. Contain your excitement.

This library seeks to implement several computational geometry tasks
in a parallel setting.

The data structure used in this library is based on the topological
edge algebras and the quad-edge structure for representing arbitrary
subdivisions of two-dimensional manifolds.

## Todo Items

Compare memory efficiency vs. problem size to theoretical bound using
valgrind.

Calculate the memory efficiency of the Delaunay class (using
quad-edge) vs a standard conforming finite element mesh
representation.
