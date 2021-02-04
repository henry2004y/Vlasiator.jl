```@meta
CurrentModule = Vlasiator
```

# Vlasiator

Data processing and analyzing tool for the numerical model for collisionless ion-kinetic plasma physics [Vlasiator](https://github.com/fmihpc/vlasiator).
This package is built upon its sister in Python [Analysator](https://github.com/fmihpc/analysator), and is carefully designed for performance and ease of use.

It contains the following features:
* Reading [VLSV](https://github.com/fmihpc/vlsv) format data.
* Extracting quantities from the simulation at a given point/line/cut.
* Plotting 2D cuts from 2D/3D simulations.

!!! note
    This package is still young, so be careful for any future breaking changes!

## Getting started

To install it,
```
pkg> add Vlasiator
```

You can then get started with
```
julia> using Vlasiator
```

More usages can be found in the [examples](examples.md).