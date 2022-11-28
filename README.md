# Vlasiator.jl

<p align="center">
  <img src="docs/src/figures/logo_fancy_black.png" height="200"><br>
  <a href="https://github.com/henry2004y/Vlasiator.jl/actions">
    <img src="https://img.shields.io/github/workflow/status/henry2004y/Vlasiator.jl/CI">
  </a>
  <a href="https://codecov.io/gh/henry2004y/Vlasiator.jl">
    <img src="https://img.shields.io/codecov/c/github/henry2004y/Vlasiator.jl">
  </a>
  <a href="https://henry2004y.github.io/Vlasiator.jl/stable">
    <img src="https://img.shields.io/badge/docs-stable-blue">
  </a>
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue">
  </a>
</p>

Data processing and analyzing tool for the collisionless ion-kinetic plasma physics numerical model [Vlasiator](https://github.com/fmihpc/vlasiator).

## Installation

In the Julia REPL,

```julia
julia> ]
pkg> add Vlasiator
```

Visualization via [PyPlot](https://github.com/JuliaPy/PyPlot.jl), [Makie.jl](https://makie.juliaplots.org/stable/), and [Plots.jl](https://docs.juliaplots.org/stable/) are supported. Please refer to [the manual](https://henry2004y.github.io/Vlasiator.jl/stable/#Getting-started) for installing different plotting backends.

## Usage

First import the package

```julia
julia> using Vlasiator
```

For Vlasiator data, e.g. `demo.vlsv`, we can load via

```julia
julia> meta = load("demo.vlsv")
```

For plotting 2D contours with PyPlot,

```julia
julia> pcolormesh(meta, "proton/vg_rho")
```

More usages can be found in [the manual](https://henry2004y.github.io/Vlasiator.jl/stable/manual/).
