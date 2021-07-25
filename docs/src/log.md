# Log

## Precision

For post-processing and data analysis purposes, it makes less sense to stick to double precisions, so we consistently use `Float32` in Vlasiator.jl. (Not exactly, in fact in VDF plots we are still using double precision if `f` is saved in `Float64`.)

## Performance

The VLSV loader inherits the basic structure from [Analysator](https://github.com/fmihpc/analysator) and is redesigned for performance.

* Besides the language difference in speed, one of the key decisions in boosting performance is to avoid the usage of dictionary with integer keys as much as possible.
* It is generally faster to read a bunch of cell IDs together than to read each cell one-by-one.

For development, it is recommended to use [PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl) to run the test suite:
```
using PkgBenchmark, Vlasiator
results = benchmarkpkg(Vlasiator)
```
or if you want to compare the current status of the package against a different git version
```
judge(Vlasiator, "97e3dca6b2474d7bdc5b62b5bf98ecf070516e5e")
```
To export results to markdown format,
```
export_markdown("testresult", results)
```
See more in the PkgBenchmark [manual](https://juliaci.github.io/PkgBenchmark.jl/dev/).

### Benchmarks

* Reading DCCRG grid variables
| Julia | tmean [μs] | Python | tmean [μs] |
|:------|:----------:|:-------|:----------:|
| 2MB  | 200 | 2MB  | 1000   |
| 50MB | 400 | 50MB | 1000   |

* Reading field solver grid[^1] variables
| 26GB   | tmean [s] |
|:-------|:---------:|
| Julia  | 13   |
| Python | 45   |

[^1]: The field solver grid is a regular Cartesian grid at the finest refinement level. Therefore the storage requirement for fsgrid variables are quite significant: with 16 GB memory it is barely enough to read `fg_b` once; it will go out of memory for the second time! This reading time corresponds to 35% of the maximum sequential read speed on the target machine.

* From starting Julia/Python to the first plot[^2]
| 2.3GB  | tmean [s] |
|:-------|:---------:|
| Julia  | 11.6  |
| Python | 9.3   |

[^2]: This inefficieny is a famous problem in Julia known as "time to first plot". On the Python side, however, I don't know why using Analysator is slower (2.3GB file, 4.8s) than directly calling matplotlib functions (2.3GB file, 0.5s).

* Reading and plotting one 2d slice of proton density out of 3D AMR data

| 26GB   | tmean [s] |
|:-------|:---------:|
| Julia  | 0.35  |
| Python | 1.7   |

* Virtual satellite tracking from 845 frames of 3D AMR data (26G per frame) on a cluster

| 1 CPU   | tmean [m][^3] |
|:-------|:---------:|
| Julia  | 11    |
| Python | 125   |

[^3]: The timings are for a single CPU on Vorna, a local cluster at University of Helsinki with Intel Xeon CPUs. With multithreading, the Julia timings can scale linearly on a node with the number of cores used. For example, with 8 threads, Julia takes ~80s to finish.