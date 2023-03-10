# ---
# title: Virtual satellite
# id: demo_virtual_sat
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to extract data for a virtual satellite
# ---

This demo shows how to extract data for a virtual satellite. Outputs are stored in CSV format for sharing data.

Usage:
```shell
julia -t nthreads demo_virtual_satellite.jl
```
or
```shell
JULIA_NUM_THREADS=nthreads julia demo_virtual_satellite.jl
```

```julia
using Glob, DelimitedFiles, Vlasiator, DataFrames
using Vlasiator: RE # Earth radius [m]

function extract_vars(files, loc)
   nfiles = length(files)
   # variables to be extracted
   t   = zeros(Float32, nfiles)
   rho = zeros(Float32, nfiles)
   v   = zeros(Float32, 3, nfiles)
   p   = zeros(Float32, nfiles)
   b   = zeros(Float32, 3, nfiles)
   e   = zeros(Float32, 3, nfiles)

   id = load(files[1]) do meta
      getcell(meta, loc)
   end

   # Extract data from each frame
   Threads.@threads for i = eachindex(files)
      meta = load(files[i])
      t[i] = meta.time
      rho[i] = readvariable(meta, "proton/vg_rho", id)[1]
      v[:,i] = readvariable(meta, "proton/vg_v", id)
      p[i] = readvariable(meta, "vg_pressure", id)[1]
      b[:,i] = readvariable(meta, "vg_b_vol", id)
      e[:,i] = readvariable(meta, "vg_e_vol", id)
   end

   df = DataFrame(t = t, rho = rho, v = v, p = p, b = b, e = e)
   # Save into text file
   writedlm("satellite.csv", Iterators.flatten(([names(df)], eachrow(df))), ',')
end

#####

files = glob("bulk*.vlsv", "./")

# virtual satellite location
const loc = Float64[12, 0, 0] .* RE

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
```