# Log

## Test Data

If you don't have VLSV data at hand, Vlasiator.jl provides some test data for you to begin with.

```
using LazyArtifacts

rootpath = artifact"testdata"
files = joinpath.(rootpath, ("bulk.1d.vlsv", "bulk.2d.vlsv", "bulk.amr.vlsv"))
```

These are also used in the standard test. These will be automatically downloaded from [vlsv_data](https://github.com/henry2004y/vlsv_data) if you run the package test locally.

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

!!! note
    The numbers shown here are comparisons between Analysator v0.9 and Vlasiator.jl v0.9.22 running Python 3.6.9 and Julia 1.7.3. The timings are performed on a i5-10210U @ 1.6GHz if not specified. Keep in mind that when we are comparing against Python, we are mostly likely comparing with the underlying C libraries with a Python wrapper.

* Reading DCCRG grid variables
| Variable | 80KB Float32 | 900KB Float64 | 32MB Float64 |
|:------|:----------:|:-------|:----------:|
| Julia  [ms] | 0.02 | 0.4 | 12.5[^1] |
| Python [ms] | 0.14 | 8.7  | 295 |

[^1]: Julia can be even faster if there is no conversion from Float64 to Float32. See [Precision](#precision).

* Reading field solver grid variables[^2]
| 13 GB  | tmean [s] |
|:-------|:---------:|
| Julia  | 8   |
| Python | 61  |

[^2]: The field solver grid is a regular Cartesian grid at the finest refinement level. Therefore the storage requirement for fsgrid variables are quite significant: with 16 GB memory it is barely enough to read `fg_b` once. It will go out of memory for the second time in Analysator, but not in Vlasiator.jl --- see [Memory](#memory). This reading time corresponds to 35% of the maximum sequential read speed on the target machine.

* From starting Julia/Python to the first plot of 2D density contour[^3]
| 28 MB  | tmean [s] |
|:-------|:---------:|
| Julia 1.7  | 11.0  |
| Python 3.6 | 9.5   |

[^3]: This inefficieny of Julia is a famous problem in the community known as "time to first plot". On the Python side, however, I don't know why using Analysator is slower (2.3GB file, 4.8s) than directly calling matplotlib functions (2.3GB file, 0.5s).

* Reading and plotting one 2d slice of proton density out of 3D AMR data

| 32 MB  | tmean [s] |
|:-------|:---------:|
| Julia  | 0.35  |
| Python | 1.7   |

* Static virtual satellite tracking from 200 frames of 3D AMR data (26G per frame, 32 MB Cell IDs) on a cluster

| 1 core | tmean [s][^4] |
|:-------|:---------:|
| Julia  | 173   |
| Python | 316   |

[^4]: A single CPU core is used on Vorna, a local cluster at University of Helsinki with Intel Xeon E5-2697 @ 2.70GHz. With multithreading, the Julia timings can scale linearly on a node with the number of cores used.

## Precision

For post-processing and data analysis purposes, it makes less sense to stick to double precisions, so we mostly use `Float32` in Vlasiator.jl. Several exceptions are:

* physical constants are defined in `Float64`, since single precision only resolve up to ±3.4E+38, and it may go out of bound in the middle of calculation (e.g. plasma frequency).

## Int v.s. UInt

We have not made a consensus on which integer to use for cell indexes. Be careful about potential bugs due to incorrect arithmetics especially for unsigned integers!

## Memory

Vlasiator output files can be large. If we have limited memory relative to the file size, Vlasiator.jl provide direct hard disk mapping through `mmap` in Julia. With this mechanism you never need to worry about unable to process data with small free memory.

## Parallelism

The current design choice is to achieve optimal serial performance per file, and apply parallel processing across individual files. In most common cases, the time it takes for post-processing one snapshot is reasonably short, but the number of snapshots are large. Julia's built-in support for all kinds of parallelism paradigm (multithreads, multiprocess, channel) and external support from packages (MPI.jl, Polyester.jl) can be relatively easily incorported to make the whole workflow parallel.

In the [examples](https://github.com/henry2004y/Vlasiator.jl/tree/master/examples), you can find the usages of

* multi-threading with `@threads` (recommended when working within one node)
* multi-processing with `pmap` 
* multi-processing with `RemoteChannel`
* `ClusterManagers` for multi-node jobs

## VTK

VLSV is just an uncompressed binary format. If we convert VLSV to VTK through `write_vtk`, the generated VTK files, even the highest resolution one with every coarse cell mapping to the finest level, can be several times smaller than the original VLSV file.

One drawback of this conversion is that it cannot deal with phase space outputs, i.e. VDFs.