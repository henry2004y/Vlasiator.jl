# Benchmarks

The test file information are listed below:

| Index  | Filename        | Number of Cells | Dimension | AMR | Public |
|:-------|:----------------|:---------------:|:---------:|:---:|:------:|
| 1      | 1d_single.vlsv  | 20              | 1         | No  | Yes    |
| 2      | bulk.2d.vlsv    | 6,300           | 1         | No  | Yes    |
| 3      | 2d_double.vlsv  | 51,200          | 2         | No  | Yes    |
| 4      | 2d_AFC.vlsv     | 4,612,500       | 2         | No  | Yes    |
| 5      | 3d_EGI.vlsv     | 3,966,580       | 3         | Yes | No     |

Access to the public data can be found from [vlsv_data](https://github.com/henry2004y/vlsv_data).

!!! note
    The numbers shown here are comparisons between Analysator v0.9 and Vlasiator.jl v0.9.32 running Python 3.6.9/3.9.7 and Julia 1.8.3 with the scripts [`perf.jl`](https://github.com/henry2004y/Vlasiator.jl/blob/master/benchmark/perf.jl) and [`perf.py`](https://github.com/henry2004y/Vlasiator.jl/blob/master/benchmark/perf.py). The timings are collected from a i5-10210U @ 1.6GHz CPU with 16 GB RAM if not specified.

* Loading meta data[^1]
| File Index | Julia [ms] | Python [ms] | Ratio |
|:-----------|:----------:|:-----------:|:-----:|
| 1          | 0.18       | 1.19        | 6.6   |
| 2          | 0.51       | 1.66        | 3.1   |
| 3          | 2.33       | 3.11        | 1.3   |
| 4          | 506        | 277         | 0.5   |
| 5          | 549        | 283         | 0.5   |

[^1]: See the [issue](https://github.com/henry2004y/Vlasiator.jl/issues/124) about sorting. The performance of EzXML is also not ideal: we may need to find a better XML parser in Julia.

* Reading DCCRG grid variables
| File Index | Julia [ms] | Python [ms] | Ratio  |
|:-----------|:----------:|:-----------:|:------:|
| 1          | 0.004      | 0.07        | 17     |
| 2          | 0.02[^2]   | 0.08        | 4      |
| 3          | 0.17[^2]   | 0.21        | 1.2    |
| 4          | 20[^2]     | 23          | 1.1    |
| 5          | 11[^2]     | 11          | 1.0    |

[^2]: Vlasiator.jl can be faster if there is no conversion from Float64 to Float32. See [Precision](log.md#precision).

* Reading field solver grid variables[^3]
| File Index | Size            | Julia [s] | Julia, mmap [s] | Python [s] | Ratio |
|:-----------|:----------------|:---------:|:---------------:|:----------:|:-----:|
| 5          | 6.2 GiB Float32 | 18        | 8               | 56         | 7     |

[^3]: The field solver grid is a regular Cartesian grid at the finest refinement level introduced after Vlasiator 5. Therefore fsgrid variables are quite large for 3D AMR runs: with limited memory (e.g. 16 GB RAM) you may encounter out-of-memory issues when reading `fg_b` more than once. In Vlasiator.jl, we provide the option `usemmap=true` for reading large arrays --- see [Memory](log.md#memory) for more.

* Plotting 2D density contours on a uniform mesh (no GUI)
| File Index | Julia [s] | Python [s] | Ratio |
|:-----------|:---------:|:----------:|:-----:|
| 4          | 0.5       | 5.4       | 11   |

* Plotting 2D density slices from a 3D AMR mesh (no GUI)
| File Index | Julia [s][^4] | Python [s] | Ratio |
|:-----------|:-------------:|:----------:|:-----:|
| 5          | 0.5           | 3.1        | 6.2    |

[^4]: The first time execution will be slower due to JIT compilation (which is excluded in the timing here). This is known as "Time-To-First-X" in the Julia community.

* Static virtual satellite tracking from 3D AMR data (26G per frame, 32 MB Cell IDs) on a cluster[^5]

| Frames | Julia, 1 thread [s] | Julia, 2 threads [s] | Python [s] |
|:-------|:-------------------:|:--------------------:|:----------:|
| 845    | 45                  | 34                   | 1220       |

[^5]: University of Helsinki cluster Vorna with Intel Xeon E5-2697 @ 2.70GHz.