```@meta
CurrentModule = Vlasiator
```

# Vlasiator.jl

Data processing and analyzing tool for the numerical model for collisionless ion-kinetic plasma physics [Vlasiator](https://github.com/fmihpc/vlasiator).
This lightweight package is built upon its sister in Python [Analysator](https://github.com/fmihpc/analysator) and carefully designed for performance, capability and ease of use.
It can be easily integrated with external packages like [FieldTracer.jl](https://github.com/henry2004y/FieldTracer.jl) and [TestParticle.jl](https://github.com/henry2004y/TestParticle.jl) to do all kinds of in-depth analysis.

Vlasiator.jl contains the following features:

* Reading [VLSV](https://github.com/fmihpc/vlsv) format data.
* Calculating derived quantities from VLSV outputs.
* Extracting quantities at a given point/line/plane.
* Plotting 1D curves/2D cuts of saved/derived variables, and phase space distributions.
* Converting VLSV into VTK format for postprocessing in e.g. ParaView and VisIt.

!!! warning
    This package mostly aims at supporting Vlasiator 5.0+. Older versions of Vlasiator has different naming standard for outputs, and is not guaranteed to work. [This analysator wiki page](https://github.com/fmihpc/analysator/wiki/Supported-variables-and-data-reducers) describes the old and new naming standards in detail.

## Getting started

To install,

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