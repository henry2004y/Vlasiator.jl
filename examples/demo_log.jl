# Example of extracting log file timing information and visualize them.

using Vlasiator, Dates, Plots

# plotly is nice for scanning through data interactively
plotly()

const file = "logfile.txt"

timestamps, speed = readlog(file)

scatter(timestamps, speed,
   markershape=:circle,
   #yaxis=:log10, 
   xlabel="Time", ylabel="Time per simulated second [s]")