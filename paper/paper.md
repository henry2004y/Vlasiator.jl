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
date: 14 Semptember 2021
bibliography: paper.bib

---

# Summary

Space weather is used to describe the environmental effects in the solar system
caused by the solar wind, a stream of charged particles carrying the solar
electromagnetic field. Vast majority of space in the solar system is filled with
charged particles, i.e. plasma.  Plasma can carry electromagnetic field and
interacts with astronomical object's magnetic field to create a magnetosphere
near the object. `Vlasiator` [@palmroth2018vlasov] is a numerical model for
collisionless ion-kinetic plasma physics, aiming at studying space weather in
the global magnetosphere.

Due to the multi-dimensional approach at ion scales, Vlasiator's computational
challenges are immense. The storage required to resolve the phase space
distributions can easily go beyond tegabytes with each reduced snapshot goes
beyond 10 GB, which requires efficient numerical tools for processing the data.

# Statement of need

`Vlasiator.jl` is a Julia package for processing and analyzing simulation data
from Vlasiator.
This lightweight package is built upon its sister package in Python `Analysator`
[@analysator] and is carefully designed for performance, capability and ease of
use. It can be easily integrated with other open source packages in the
community like `FieldTracer.jl` [@fieldtracer] for tracing along the field lines
and `TestParticle.jl` [@testparticle] for test particle simulations.

`Vlasiator.jl` contains the following main features:

- Reading `VLSV` format data [@vlsv], including `DCCRG` [@honkonen2013parallel]
and `FSGRID` [@fsgrid], at any size.
- Calculating derived quantities from VLSV outputs.
- Extracting quantities at a given point/line/plane.
- Plotting 1D curves/2D cuts of saved variables, and phase space distributions.
- Converting VLSV into VTK format for postprocessing in e.g. ParaView and VisIt.

`Vlasiator.jl` is targeted at space physics researchers who want to visualize
and analyze Vlasiator simulation outputs in an efficient manner.
It achieves optimal serial performance for single file processing and can be
directly applied to parallel batch jobs using both multithreads and
multiprocesses. It has preliminarily been used in ultra-low frequency wave
studies under time-varying solar wind conditions [@ressac].

The combination of speed and design in `Vlasiator.jl` will enable exciting
scientific explorations of forthcoming data from the large scale simulation by
students and experts alike.

# Acknowledgements

We acknowledge support from Lucile Turc during the genesis of this project.
Funding for this work is provided by Finnish Institute of Science.

# References