---
title: 'Vlasiator.jl: A Julia package for processing Vlasiator data'
tags:
  - Julia
  - space physics
  - plasma
authors:
  - name: Hongyang Zhou
    orcid: 0000-0003-4571-4501
    affiliation: "1"
affiliations:
 - name: University of Helsinki
   index: 1
date: 28 November 2022
bibliography: paper.bib
---

# Summary

`Vlasiator.jl` is a Julia [@julia] package for processing and analyzing simulation data from the collisionless ion-kinetic plasma physics numerical model [`Vlasiator`](https://github.com/fmihpc/vlasiator) [@vlasiator5.2.1].
This lightweight package retains all the actively used funtionalities in its sister Python package [`Analysator`](https://github.com/fmihpc/analysator) [@analysator] and is carefully designed for performance, capability and ease of use.

`Vlasiator.jl` contains the following main features:

- Reading [VLSV](https://github.com/fmihpc/vlsv) format data, including [`DCCRG`](https://github.com/fmihpc/dccrg) [@honkonen2013parallel] and [FSGRID](https://github.com/fmihpc/fsgrid), at any size.
- Calculating predefined derived quantities from raw VLSV outputs.
- Extracting quantities at a given point/line/plane/box.
- Visualizing 1D curves/2D cuts/3D volumes of saved/derived variables and phase space distributions via multiple visualization libraries such as `Matplotlib` [@matplotlib], `Plots.jl` [@plots], and Makie.jl [@makie].
- Analyzing the velocity distribution functions reconstructed from sparsity storage.
- Converting the selected part or whole data from VLSV into VTK format for post-processing and 3D rendering in ParaView and VisIt.

`Vlasiator.jl` achieves optimal serial performance for single file processing and can be directly applied to parallel batch jobs using both multithreads and multiprocesses.
The interoperability with Python can be easily achieved via two community packages `JuliaCall` [@juliacall] and `PyJulia` [@pyjulia].
It can also be easily integrated with other packages in the community like [FieldTracer.jl](https://github.com/henry2004y/FieldTracer.jl) for tracing along the field lines and [TestParticle.jl](https://github.com/henry2004y/TestParticle.jl) for embedded test particle simulations.
The performance and ease-of-use of `Vlasiator.jl` will enable exciting reproducible scientific exploration of forthcoming data from exascale simulations.

# Statement of need

Space weather is used to describe the environmental effects in the solar system caused by the solar wind, a stream of charged particles carrying the solar electromagnetic field.
Vast majority of space in the solar system is filled with charged particles, i.e. plasma.  Plasma can carry electromagnetic field and interacts with astronomical object's magnetic field to create a magnetosphere near the object.
`Vlasiator` [@palmroth2018vlasov] is a numerical model for collisionless ion-kinetic plasma physics, aiming at studying space weather in the global magnetosphere.
Due to the multi-dimensional approach at ion scales, `Vlasiator`'s computational challenges are immense.
The storage required to resolve the phase space distributions can easily go beyond terabytes with each reduced snapshot goes beyond ~10 GB, which necessitates the development of high performance programs for processing the data.

`Vlasiator.jl` tackles the post-processing challenges by taking advantage of novel techniques shared in the open source community, which is built from the ground up to leverage the power of Julia and successful tools written in C++ and Python.
It is targeted at space plasma physics researchers who want to analyze and visualize `Vlasiator` simulation outputs in an efficient manner: we have reached 8 - 55 times speedups over the same tasks in `Analysator`.
This package satisfies the current requirements of simulation data processing, unifies the implementation in a single language base, and facilitates more fluent downstream processing tasks and in-depth exploration of large datasets.
It has been used extensively in ultra-low frequency wave studies [@pc5] and responses of near-Earth space under changing solar wind conditions [@ressac].

# Acknowledgements

We acknowledge the support from Lucile Turc during the genesis of this project.
Funding for this work is provided by the University of Helsinki (three-year research grant 2020-2022) and the Academy of Finland (grant numbers 322544 and 328893).

# References