# ---
# title: Log tracking
# id: demo_log
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to extract logfile.txt information
# ---

This is an example of extracting log file timing information and visualize them.
```julia
using Vlasiator, Dates, Plots

# plotly is nice for scanning through data interactively
plotly()

const file = "logfile.txt"

timestamps, speed = readlog(file)

scatter(timestamps, speed,
   markershape=:circle,
   #yaxis=:log10,
   xlabel="Time", ylabel="Time per simulated second [s]")
```