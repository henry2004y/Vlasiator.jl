# FAQ

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
