```@meta
CurrentModule = Vlasiator
```

# Vlasiator.jl

Data processing and analyzing tool for the numerical model for collisionless ion-kinetic plasma physics [Vlasiator](https://github.com/fmihpc/vlasiator).
This lightweight package is built upon its sister in Python [Analysator](https://github.com/fmihpc/analysator), and is carefully designed for performance and ease of use.
It becomes even more powerful when using together with other fantastic packages: combined with external packages like [FieldTracer.jl](https://github.com/henry2004y/FieldTracer.jl) and [TestParticle.jl](https://github.com/henry2004y/TestParticle.jl), it is possible to do all kinds of in-depth analysis.

Vlasiator.jl contains the following features:
* Reading [VLSV](https://github.com/fmihpc/vlsv) format data.
* Converting VLSV into VTK format.
* Extracting quantities from the simulation at a given point/line/cut.
* Plotting 2D cuts.
* Plotting phase space distributions.

!!! note
    This package is still young, so be careful for any future breaking changes!

!!! warning
    This package mostly aims at supporting Vlasiator 5.0+. Older versions of Vlasiator has different naming standard for outputs, and is not guaranteed to work.

## Getting started

To install it,
```
pkg> add Vlasiator
```

You can then get started with
```
julia> using Vlasiator
```

If you want to use Plots.jl for visualization, add it also through the pkg manager; if you aim at using Matplotlib, besides adding `PyPlot`, you should also link to a preinstalled Python version by setting the environment variable and building the PyCall package
```
ENV["PYTHON"]="your python executable"
Pkg.build("PyCall")
```
Details are described in [automated matplotlib installation](https://github.com/JuliaPy/PyPlot.jl#automated-matplotlib-installation).

## Author

This module is written by Hongyang Zhou.