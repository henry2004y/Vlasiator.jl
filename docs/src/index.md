```@meta
CurrentModule = Vlasiator
```

# Vlasiator.jl

Data processing and analyzing tool for the numerical model for collisionless ion-kinetic plasma physics [Vlasiator](https://github.com/fmihpc/vlasiator).
This lightweight package is built upon its sister in Python [Analysator](https://github.com/fmihpc/analysator) and carefully designed for performance, capability and ease of use.
It can be easily integrated with external packages like [FieldTracer.jl](https://github.com/henry2004y/FieldTracer.jl) and [TestParticle.jl](https://github.com/henry2004y/TestParticle.jl) to do all kinds of in-depth analysis.

Vlasiator.jl contains the following features:

- Reading [VLSV](https://github.com/fmihpc/vlsv) format data.
- Calculating derived quantities from raw VLSV outputs.
- Extracting quantities at a given point/line/plane/box.
- Plotting 1D curves/2D cuts of saved/derived variables and phase space distributions.
- Analyzing velocity distribution functions.
- Appending [DCCRG](https://github.com/fmihpc/dccrg) arrays to VLSV files.
- Converting selected domain and variables from VLSV into VTK format for data analysis and visualization in ParaView and VisIt.
- Monitoring Vlasiator run log files.

!!! warning
    This package mostly aims at supporting Vlasiator 5.0+. Older versions of Vlasiator has different naming standard for outputs, and is not guaranteed to work. [This analysator wiki page](https://github.com/fmihpc/analysator/wiki/Supported-variables-and-data-reducers) describes the old and new naming standards in detail.

## Getting started

To install,

```julia
julia> ]
pkg> add Vlasiator
```

You can then use the package via

```julia
julia> using Vlasiator
```

If you want to use [Plots.jl](https://docs.juliaplots.org/stable/) or [Makie.jl](https://makie.juliaplots.org/stable/) for visualization, add them through the pkg manager; if you aim at using Matplotlib, besides adding [`PyPlot`](https://github.com/JuliaPy/PyPlot.jl), you should also link to a preinstalled Python version by setting the environment variable and building the PyCall package

```julia
ENV["PYTHON"]="your python executable"
Pkg.build("PyCall")
```

If `ENV["PYTHON"] = ""` before building, a private Python distribution will be installed via Miniconda. Details are described in [automated matplotlib installation](https://github.com/JuliaPy/PyPlot.jl#automated-matplotlib-installation).

## Author

This module is written by Hongyang Zhou.