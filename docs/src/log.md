# Log

## Test Data

If you don't have VLSV data at hand, Vlasiator.jl provides some test data for you to begin with.

```julia
using LazyArtifacts

rootpath = artifact"testdata"
files = joinpath.(rootpath, ("bulk.1d.vlsv", "bulk.2d.vlsv", "bulk.amr.vlsv"))
```

These are also used in the standard test. These will be automatically downloaded from [vlsv_data](https://github.com/henry2004y/vlsv_data) if you run the package test locally.

## Performance

The VLSV loader inherits the basic structure from [Analysator](https://github.com/fmihpc/analysator) and is redesigned for performance.

* For general data reading, a dictionary is constructed for cell IDs and orderings for O(1) timings.
* It is faster to read a bunch of cell IDs together, if possible, than to read each cell one-by-one.
* Specific methods are provided for targeted tasks that are much faster than the generic approaches. For instance, `extractsat` for multi-frame static satellite extraction can be more than 10x faster than first reading the metadata and then extracting variables for each frame.

For development, it is recommended to use [PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl) to run the test suite:

```julia
using PkgBenchmark, Vlasiator
results = benchmarkpkg(Vlasiator)
```

or if you want to compare the current status of the package against a different git version

```julia
judge(Vlasiator, "97e3dca6b2474d7bdc5b62b5bf98ecf070516e5e")
```

To export results to Markdown,

```julia
export_markdown("testresult", results)
```

See more in the PkgBenchmark [manual](https://juliaci.github.io/PkgBenchmark.jl/dev/).

### Benchmarks

!!! note
    The numbers shown here are comparisons between Analysator v0.9 and Vlasiator.jl v1.0 running Python 3.6.9/3.9.7 and Julia 1.8.0-rc3. The timings are performed on a i5-10210U @ 1.6GHz if not specified. Note that when we are comparing against Python, we are mostly likely comparing with the underlying C libraries with a Python wrapper.

* Reading DCCRG grid variables
| Variable | 80KB Float32 | 900KB Float64 | 32MB Float64 |
|:------|:----------:|:-------|:----------:|
| Julia 1.8  [ms] | 0.02 | 0.4 | 13[^1] |
| Python 3.6 [ms] | 0.14 | 8.7 | 295 |

[^1]: Vlasiator.jl can be even faster if there is no conversion from Float64 to Float32. See [Precision](#precision).

* Reading field solver grid variables[^2]
| 13 GB  | tmean [s] |
|:-------|:---------:|
| Julia 1.8  | 8   |
| Python 3.6 | 61  |

[^2]: The field solver grid is a regular Cartesian grid at the finest refinement level. Therefore the storage requirement for fsgrid variables are quite significant: with 16 GB memory it is barely enough to read `fg_b` once. It will go out of memory for the second time in Analysator, but not in Vlasiator.jl --- see [Memory](#memory). This reading time corresponds to 35% of the maximum sequential read speed on the target machine.

* Plotting 2D density contours on a uniform mesh[^3]
| 28MB | tmean [s] |
|:-------|:---------:|
| Julia 1.8, 1st time | 2.4  |
| Julia 1.8, 2nd time | 0.5  |
| Python 3.6 | 4.7   |

* Plotting 2D density slices from an 3D AMR mesh
| 32MB | tmean [s] |
|:-------|:---------:|
| Julia 1.8, 1st time | 2.7  |
| Julia 1.8, 2nd time | 0.5  |
| Python 3.6 | 5.0  |

[^3]: The inefficieny of JIT for the first time execution is a famous problem in the Julia community known as "Time-To-First-X".

* Static virtual satellite tracking from 200 frames of 3D AMR data (26G per frame, 32 MB Cell IDs) on a cluster

| 1 core | tmean [s][^4] |
|:-------|:---------:|
| Julia 1.8  | 14    |
| Python 3.9 | 302   |

[^4]: A single CPU core is used on Vorna, a local cluster at University of Helsinki with Intel Xeon E5-2697 @ 2.70GHz. With multithreading, the Julia timings can scale linearly on a node with the number of cores used.

## Precision

For post-processing and data analysis purposes, it makes less sense to stick to double precisions, so we mostly use `Float32` in Vlasiator.jl for numerical arrays. Several exceptions are:

* physical constants are defined in `Float64`, since single precision only resolves up to ±3.4E+38, and it may go out of bound in the middle of calculation (e.g. plasma frequency).

## Int vs. UInt

Integers but not unsigned integers shall be used for indexing, even though [unsigned integers are tempting](http://eigen.tuxfamily.org/index.php?title=FAQ#Why_Eigen.27s_API_is_using_signed_integers_for_sizes.2C_indices.2C_etc..3F).

## Memory

Vlasiator output files can be large. If we have limited memory relative to the file size, Vlasiator.jl provide direct hard disk mapping through [`mmap`](https://docs.julialang.org/en/v1/stdlib/Mmap/) in Julia. With this mechanism you never need to worry about unable to process data with small free memory. Besides, we found that proper usage of `mmap` can also speed up reading and reduce memory comsumption. However, without `reinterpret` we may encounter the alignment issue!

## Parallelism

The current design choice is to achieve optimal serial performance per file, and apply parallel processing across individual files. In most common cases, the time it takes for post-processing one snapshot is reasonably short, but the number of snapshots are large. Julia's built-in support for all kinds of parallelism paradigm (multithreading, multiprocessing, channel) and external support from packages (MPI.jl, Polyester.jl) can be relatively easily incorported to make the whole workflow parallel.

* multi-threading with `@threads` (recommended when working within one node)
* multi-processing with `pmap` 
* multi-processing with `RemoteChannel`
* `ClusterManagers` for multi-node jobs

See more in the [examples](https://github.com/henry2004y/Vlasiator.jl/tree/master/examples).

## VTK

VLSV is just an uncompressed binary format. If we convert VLSV to VTK through `write_vtk`, the generated VTK files, even the highest resolution one with every coarse cell mapping to the finest level, can be several times smaller than the original VLSV file.

One drawback of this conversion is that it cannot deal with phase space outputs, i.e. VDFs.