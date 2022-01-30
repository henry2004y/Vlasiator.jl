# Sample postprocessing script for virtual satellite tracking.
# Outputs are stored in CSV format for sharing data.
#
# Usage:
#   julia -t nthreads demo_virtual_satellite.jl
# or
#   JULIA_NUM_THREADS=nthreads julia demo_virtual_satellite.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Glob, DelimitedFiles, Vlasiator, DataFrames
using Vlasiator: RE # Earth radius [m]

function extract_vars(files, loc)
   nfiles = length(files)
   # variables to be extracted
   t   = zeros(Float32, nfiles)
   rho = zeros(Float32, nfiles)
   vx  = zeros(Float32, nfiles)
   vy  = zeros(Float32, nfiles)
   p   = zeros(Float32, nfiles)
   bz  = zeros(Float32, nfiles)
   ex  = zeros(Float32, nfiles)
   ey  = zeros(Float32, nfiles)

   local id
   let meta = load(files[1])
      id = getcell(meta, loc)
   end

   # Extract data from each frame
   Threads.@threads for i = eachindex(files)
      meta = load(files[i])
      t[i] = meta.time
      rho[i] = readvariable(meta, "proton/vg_rho", id)[1]
      vx[i], vy[i] = readvariable(meta, "proton/vg_v", id)[1:2]
      p[i] = readvariable(meta, "vg_pressure", id)[1]
      bz[i] = readvariable(meta, "vg_b_vol", id)[3]
      ex[i], ey[i] = readvariable(meta, "vg_e_vol", id)[1:2]
   end

   df = DataFrame(t = t, rho = rho, vx = vx, vy = vy, p = p, bz = bz, ex = ex, ey = ey)
   # Save into text file
   writedlm("satellite.csv", Iterators.flatten(([names(df)], eachrow(df))), ',')
end

#####

files = glob("bulk*.vlsv", "./")

# virtual satellite location
loc = [12Re, 0, 0]

println("Number of files: $(length(files))")
println("Extracting location: $loc")
println("Running with $(Threads.nthreads()) threads...")

@time extract_vars(files, loc)

println("Virtual satellite extraction done!")

## Visualization
#=
using PyPlot, DelimitedFiles

data = readdlm("satellite.csv", ','; header=true)

fig, ax = subplots(figsize=(8,10), 5,1, sharex=true, constrained_layout=true)

ax[1].plot(data[1][:,1], data[1][:,2] ./ 1e6, label="density")
ax[2].plot(data[1][:,1], data[1][:,3] ./ 1e3, label="vx")
ax[2].plot(data[1][:,1], data[1][:,4] ./ 1e3, label="vy")
ax[3].plot(data[1][:,1], data[1][:,5] .* 1e9, label="p")
ax[4].plot(data[1][:,1], data[1][:,6] .* 1e9, label="bz")
ax[5].plot(data[1][:,1], data[1][:,7] .* 1e3, label="ex")
ax[5].plot(data[1][:,1], data[1][:,8] .* 1e3, label="ey")

for a in ax
   a.grid(true)
   a.legend()
end

ax[1].set_ylabel("density [amu/cc]", fontsize=14)
ax[2].set_ylabel("velocity [km/s]", fontsize=14)
ax[3].set_ylabel("pressure [nPa]", fontsize=14)
ax[4].set_ylabel("magnetic field [nT]", fontsize=14)
ax[5].set_ylabel("electric field [mV/m]", fontsize=14)
ax[5].set_xlabel("time [s]", fontsize=14)

fig.suptitle("Density Pulse Run, location = [12, 0, 0]", fontsize="xx-large")

savefig("virtual_satellite.png", bbox_inches="tight")
=#