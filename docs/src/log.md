# Log

## Performance

The VLSV loader inherits the basic structure from [Analysator](https://github.com/fmihpc/analysator) and is redesigned for performance.

* Besides the language difference in speed, one of the key decisions in boosting performance is to avoid the usage of dictionary with integer keys as much as possible.
* It is generally faster to read a bunch of cell IDs together than to read each cell one-by-one.

For development, it is recommended to use [PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl) to run the test suite:
```
using PkgBenchmark, Vlasiator
benchmarkpkg(Vlasiator)
```
or if you want to compare the current status of the package against a different git version
```
judge(Vlasiator, "97e3dca6b2474d7bdc5b62b5bf98ecf070516e5e")
```
See more in the PkgBenchmark manual.

### Benchmarks

Initial tests on reading variables from sample VLSV files: 

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