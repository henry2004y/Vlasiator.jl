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

## What is Vlasiator

Vlasiator solves the Vlasov–Maxwell system of equations for ions. The fundamental description of charged particle motion in an electromagnetic field is given by the Vlasov equation

```math
\frac{\partial f_\alpha}{\partial t} + \mathbf{v}_\alpha\frac{\partial f_\alpha}{\partial \mathbf{r}_\alpha} + \mathbf{a}_\alpha\cdot \frac{\partial f_\alpha}{\partial \mathbf{v}_\alpha} = 0
```

where ``\alpha`` denotes the particle species, ``\mathbf{r}_\alpha$ and $\mathbf{v}_\alpha`` are the spatial and velocity coordinates, ``f(\mathbf{r},\mathbf{v},t)`` is the six-dimensional phase-space density of a particle species with mass ``m`` and charge ``q``, and acceleration ``\mathbf{a}`` is given by the Lorentz force

```math
\mathbf{a}_\alpha = \frac{q_\alpha}{m_\alpha}(\mathbf{E}+\mathbf{v}_\alpha\times\mathbf{B})
```

where ``\mathbf{E}`` and ``\mathbf{B}`` are the electric and magnetic field, respectively.

The bulk parameters of the plasma, such as the ion density ``n_\alpha`` and velocity ``\mathbf{u}_\alpha``, are obtained as velocity moments of the ion velocity distribution function

```math
\begin{aligned}
n_\alpha &= \int f(\mathbf{r},\mathbf{v},t)d^3 v \\
\mathbf{u}_\alpha &= \int \mathbf{v}f(\mathbf{r},\mathbf{v},t)d^3 v
\end{aligned}
```

The magnetic field is updated using Faraday's law:

```math
\nabla \times \mathbf{E} = -\frac{\partial \mathbf{B}}{\partial t}
```

and the electric field is given by the generalized Ohm's law:

```math
\mathbf{E}_\alpha = -\mathbf{u}_\alpha \times\mathbf{B} + \frac{1}{n_\alpha e}(\nabla\times\mathbf{B})\times\mathbf{B} - \frac{1}{n_\alpha e}\nabla\cdot\overleftrightarrow{P}_e + \eta \mathbf{j}
```

The four terms on the right-handed side are the convection term, the Hall term, the electron pressure gradient term, and the resistive term, respectively.
The total current density ``\mathbf{j}`` is obtained from Ampère–Maxwell's law where the displacement current has been neglected:

```math
\nabla\times\mathbf{B} = \mu_0 \mathbf{j}
```

Finally, by determining the electron pressure tensor by using an appropriate equation of state, the evolution of the system can be followed in time. For example, let $\overleftrightarrow{P}_e = P_e \overleftrightarrow{I}$ where $P_e$ is the isotropic scalar electron pressure. In an isothermal process,

```math
P_e = n_e k_B T_e
```

where ``n_e \approx n_i`` and ``T_e`` is a constant. In an adiabatic process,

```math
P_e / n_e^\gamma = \text{const.}
```

For more details, please refer to [Vlasov methods in space physics and astrophysics](https://link.springer.com/article/10.1007/s41115-018-0003-2).

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