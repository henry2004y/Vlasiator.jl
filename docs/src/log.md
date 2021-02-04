# Log

This package was born when I was learning Vlasiator and its corresponding data structures.
The VLSV loader inherits the basic structure from [Analysator](https://github.com/fmihpc/analysator) and is redesigned for performance.

The function APIs are made trying to be consistent with Analysator.

The IOstream handle for the VLSV file requires some special attention.
In the current implementation, once the meta data is read, the file stays open until one explictly says `close(meta.fid)`.
On the Windows platform, it is not allowed to delete the file before the IO is closed.
However, this is allowed in Unix, so be careful.

## Plotting Philosophy

We should not take over what underlying plotting libraries like Matplotlib offers.
Users should be able to modify the figures as they wish even if they only know how to use the well-known plotting libraries.
Do not reinvent the wheel. For customized plotting, simply provide some sample scripts for the common tasks like zooming-in, change font sizes, add text boxes, etc..

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