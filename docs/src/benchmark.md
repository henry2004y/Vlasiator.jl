# Benchmarks

!!! note
    The numbers shown here are comparisons between Analysator v0.9 and Vlasiator.jl v0.9.32 running Python 3.6.9/3.9.7 and Julia 1.8.3 with the scripts [`perf.jl`](https://github.com/henry2004y/Vlasiator.jl/blob/master/benchmark/perf.jl) and [`perf.py`](https://github.com/henry2004y/Vlasiator.jl/blob/master/benchmark/perf.py). The timings are collected from a i5-10210U @ 1.6GHz CPU with 16 GB RAM if not specified.

* Reading DCCRG grid variables
| Size            | Julia [ms] | Python [ms] | Speedup |
|:----------------|:----------:|:-----------:|:-------:|
| 80 B Float32    | 0.002      | 0.14        | 55      |
| 25 KiB Float64  | 0.017      | 0.44        | 26      |
| 879 KiB Float64 | 0.3        | 8.7         | 29      |
| 30 MiB Float64  | 10.8       | 295         | 27      |

[^1]: Vlasiator.jl can be even faster if there is no conversion from Float64 to Float32. See [Precision](log.md#precision).

* Reading field solver grid variables[^2]
| Size            | Julia [s] | Python [s] | Speedup |
|:----------------|:---------:|:----------:|:-------:|
| 6.2 GiB Float32 | 8         | 61         | 7.6     |

[^2]: The field solver grid is a regular Cartesian grid at the finest refinement level. Therefore the storage requirement for fsgrid variables are quite significant: with limited memory (e.g. 16 GB RAM) you may encounter out-of-memory issues when reading `fg_b` more than once. In Vlasiator.jl, we provide the option `usemmap=true` for reading large arrays --- see [Memory](log.md#memory) for more.

* Plotting 2D density contours on a uniform mesh[^3]
| Size     | Julia, 1st time [s] | Julia, 2nd time [s] | Python [s] | Speedup |
|:---------|:-------------------:|:-------------------:|:----------:|:-------:|
| 26.7 MiB | 2.4                 | 0.5                 | 4.7        | 9.4     |

* Plotting 2D density slices from a 3D AMR mesh
| Size     | Julia, 1st time [s] | Julia, 2nd time [s] | Python [s] | Speedup |
|:---------|:-------------------:|:-------------------:|:----------:|:-------:|
| 30.5 MiB | 2.7                 | 0.5                 | 5.0        | 10      |

[^3]: The inefficieny of JIT for the first time execution is a famous problem in the Julia community known as "Time-To-First-X".

* Static virtual satellite tracking from 3D AMR data (26G per frame, 32 MB Cell IDs) on a cluster[^4]

| Frames | Julia, 1 thread [s] | Julia, 2 threads [s] | Python [s] |
|:-------|:-------------------:|:--------------------:|:----------:|
| 845    | 45                  | 34                   | 1220       |

[^4]: University of Helsinki cluster Vorna with Intel Xeon E5-2697 @ 2.70GHz.