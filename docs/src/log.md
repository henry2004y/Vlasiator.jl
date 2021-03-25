# Log

This package was born when I was learning Vlasiator and its corresponding data structures.
The VLSV loader inherits the basic structure from [Analysator](https://github.com/fmihpc/analysator) and is redesigned for performance.
One of the key decision in boosting performance is to avoid the usage of dictionary with integer keys as much as possible.

The function APIs are made trying to be consistent with Analysator.

The IOstream handle for the VLSV file requires some special attention.
In the current implementation, once the meta data is read, the file stays open until one explictly says `close(meta.fid)`.
On the Windows platform, it is not allowed to delete the file before the IO is closed.
However, this is allowed in Unix, so be careful.

Vlasiator has come up with many different names for the same quantities, so it is really hard to collect them correctly in post-processing. Therefore for the derived quantities, it is not guaranteed to work properly. The user should be able to compute more complicated quantities given the basic outputs from the VLSV file.

## Plotting Philosophy

We should not take over what underlying plotting libraries like Matplotlib offers.
Users should be able to modify the figures as they wish even if they only know how to use the well-known plotting libraries.
Do not reinvent the wheel. For customized plotting, simply provide some sample scripts for the common tasks like zooming-in, change font sizes, add text boxes, etc..
The original plotting APIs in Matplotlib are already complicated enough: instead of building a wrapper on top of them, it is considered a better approach to provide *recipes* such that the plotting package can understand what to do with the user defined types. Therefore the user can rely solely on the documentation of the plotting package to generate plots.

## Benchmarks

Initial tests on reading variables from a VLSV file: 

* DCCRG grid
| 2MB   | tmean [μs] |
|:-------|:---------:|
| Julia  | 200    |
| Python | 1000   |

| 50MB   | tmean [μs] |
|:-------|:---------:|
| Julia  | 400    |
| Python | 1000   |

* Field solver grid[^1]
| 26GB   | tmean [s] |
|:-------|:---------:|
| Julia  | 13   |
| Python | 45   |

[^1]: The field solver grid is a regular Cartesian grid at the finest refinement level. Therefore the storage requirement for fsgrid variables are quite significant: with 16 GB memory it is barely enough to read `fg_b` once; it will go out of memory for the second time!

I don't know why using Analysator is slower (2.3GB file, 4.8s) than directly calling matplotlib functions (2.3GB file, 0.5s).
Same thing for Julia costs 1.0s (first time ~8s including everything).

Reading and plotting one 2d slice of proton density out of 3D AMR data:

| 26GB   | tmean [s] |
|:-------|:---------:|
| Julia  | 0.35  |
| Python | 1.7   |

Virtual satellite tracking from 845 frames of 3D AMR data (26G per frame) on Vorna:

| 1 CPU   | tmean [m] |
|:-------|:---------:|
| Julia  | 11    |
| Python | 125   |

Note that the above timings are for a single CPU. With only one command added for multithreading, the Julia timings can be improved by n where n is the number of threads. For example, with 8 threads, Julia takes ~80s to finish.

## Conditional Dependency

There are certain packages that I don't want to include as dependencies, but instead I want to compile some glue codes if they are loaded. For example, I do not want to include any plotting packages as dependencies because either they are heavy-weighted, or incompatible with one another if one wants to switch.

There is a proposal in the Pkg manager for this, but it will come in later versions.
My first workaround is to include some additional scripts for a target plotting library, but that requires the location of the scripts, which is inconvenient for users.
The Plots package provides a lightweight [RecipesBase](http://juliaplots.org/RecipesBase.jl/dev/) library for defining behaviors for customized types. In the future if Makie is used, the same idea can be applied there.
The solution I am adapting now for PyPlot is to use [Requires.jl](https://github.com/JuliaPackaging/Requires.jl).
I am not alone.

## Parallelism

At some point I may want to try multi-threading in data processing. First I need to make sure adding threads does not affect single thread performance, and then I need to identify proper places for using threads rather than abuse threads at any place.