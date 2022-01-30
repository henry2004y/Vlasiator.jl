# Sample postprocessing script for wave-like structure spatial-temporal distribution
# searching.
#
# Usage:
#   julia -t 4 demo_wave_search_mt.jl
# or
#   JULIA_NUM_THREADS=4 julia demo_wave_search_mt.jl
#
# Procedures:
# 1. Extract variables in all cells from all snapshots.
# 2. For each cell and time interval, count time-series local peaks.
# 3. Plot the peak occurrence frequencies across the whole domain as an indicator for waves.
#
# Note:
# 1. When dealing with multiple variables, it is recommended to handle one variable at a
# time through the whole process due to memory considerations. With large number of frame
# counts, the current procedure is still memory-consuming.
# 2. It assumes uniform sampling in time.
# 3. Peak-finding is threaded, but plotting is still serial.
#
# Hongyang Zhou, hyzhou@umich.edu

using Glob, Vlasiator, PyPlot, LaTeXStrings

"Extract time series variable"
function extract_var(files, ncells, varname, component=0)
   nfiles = length(files)
   var = zeros(Float32, ncells[1], ncells[2], nfiles)

   # Extract data from each frame
   if component == 0 # scalar
      Threads.@threads for i = eachindex(files)
         meta = load(files[i])
         var[:,:,i] = meta[varname]
      end
   else # vector component
      Threads.@threads for i = eachindex(files)
         meta = load(files[i])
         var[:,:,i] = meta[varname][component,:]
      end
   end

   return var
end

"Count local maxima of vector `y` with moving box length `n`."
function countpeaks(y, n; interval=1)
   maxs = Int[]
   if interval == 1
      for i in 2:length(y)-1
         if y[i-1] < y[i] > y[i+1]
            push!(maxs, i)
         end
      end
   elseif interval == 2
      for i in 3:length(y)-2
         if y[i-2] ≤ y[i-1] < y[i] > y[i+1] ≥ y[i+2]
            push!(maxs, i)
         end
      end
   elseif interval == 3
      for i in 4:length(y)-3
         if y[i-3] ≤ y[i-2] ≤ y[i-1] < y[i] > y[i+1] ≥ y[i+2] > y[i+3]
            push!(maxs, i)
         end
      end
   else
      error("interval = $interval not implemented!")
   end
   nCounts = zeros(Int, length(y)-n+1)
   nCounts[1] = count(i->(1 ≤ i ≤ n), maxs)
   for i in 1:length(y)-n
      nCounts[i+1] = nCounts[i]
      if i ∈ maxs
         nCounts[i+1] -= 1
      end
      if i+n-1 ∈ maxs
         nCounts[i+1] += 1
      end
   end
   nCounts
end

"Check wave-like occurrence frequencies within box length `n` of output interval `dt`."
function checkwaves_sma(var, dt=0.5, n::Int=size(var,3); interval=1)
   nPeaks = zeros(Int, size(var,3)-n+1, size(var,1), size(var,2))

   Threads.@threads for j in axes(var, 2)
      for i in axes(var, 1)
         var_series = @view var[i,j,:]
         nPeaks[:,i,j] = countpeaks(var_series, n; interval)
      end
   end
   nPeaks ./ (n*dt)
end

function plot_dist(files, varnames, varnames_print, components, Δt, nboxlength;
   interval=1, nplotstride=1)
   @assert nboxlength ≥ 3 && isodd(nboxlength) "Expect odd box length ≥ 3!"
   if (local nfiles = length(files)) < nboxlength
      @warn "Set moving box length to the number of files..."
      nboxlength = nfiles
   end
   local x, y, tStart, tEnd, ncells
   let Re = Vlasiator.Re
      meta = load(files[1])
      tStart = meta.time
      ncells = meta.ncells
      x = LinRange{Float32}(meta.coordmin[1]/Re, meta.coordmax[1]/Re, meta.ncells[1])
      y = LinRange{Float32}(meta.coordmin[2]/Re, meta.coordmax[2]/Re, meta.ncells[2])
      meta = load(files[end])
      tEnd = meta.time
   end

   fig, ax = plt.subplots(figsize=(16,9))
   fontsize = 14
   vmin, vmax = 0.0, 0.5
   levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(vmin, vmax)
   norm = matplotlib.colors.BoundaryNorm(levels, ncolors=256, clip=true)
   ticks = range(vmin, vmax, length=11)

   fakedata = zeros(Float32, length(x), length(y))
   im = ax.pcolormesh(y, x, fakedata; norm)

   ax.set_aspect("equal")
   ax.set_xlabel(L"y [$R_E$]"; fontsize, weight="black")
   ax.set_ylabel(L"x [$R_E$]"; fontsize, weight="black")
   ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
   ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
   ax.grid(true, color="grey", linestyle="-")

   cb = fig.colorbar(im; ax, ticks, pad=0.02)
   cb.ax.set_ylabel("Frequency of local maxima occurrence, [#/s]"; fontsize)

   for i in eachindex(varnames)
      outdir = "../out/$(lowercase(varnames_print[i]))"
      !isdir(outdir) && mkdir(outdir)
      length(glob("spatial*.png", outdir)) == length(files) - nboxlength + 1 && continue

      # Obtain time series data
      var = extract_var(files, ncells, varnames[i], components[i])
      # Count local peak occuring frequencies at each location
      fPeaks = checkwaves_sma(var, Δt, nboxlength; interval)

      for it in 1:nplotstride:size(fPeaks,1) # Iterate over time
         outname = joinpath(outdir,
            "spatial_perturbation_distribution_$(lpad(it, 4, '0')).png")
         isfile(outname) && continue
         # Update plot
         im.set_array(fPeaks[it,:,:])
         ax.set_title("$(varnames_print[i]) Perturbation Detection, "*
            "t = $(round(tStart+(it-1)*Δt, digits=1)) ~ "*
            "$(round(tStart+(it+nboxlength-1)*Δt, digits=1))s";
            fontsize, fontweight="bold")

         savefig(outname, bbox_inches="tight")
      end
   end
end

##### Main

varnames = ["proton/vg_rho", "vg_pressure", "proton/vg_v", "proton/vg_v", "vg_b_vol",
   "vg_e_vol", "vg_e_vol"]
varnames_print = ["Density", "Thermal Pressure", "Vx", "Vy", "Bz", "Ex", "Ey"]
components = [0, 0, 1, 2, 3, 1, 2] # 0: scalar; 1: x, 2: y, 3: z
Δt = 0.5                           # output time interval [s]
nboxlength = 101                   # moving box average length
interval = 2                       # local peak gap minimal interval
nplotstride = 50                   # plot intervals in frames 
dir = "./" # data directory

files = glob("bulk*.vlsv", dir)

println("Total number of snapshots: $(length(files))")
println("Running with $(Threads.nthreads()) threads...")

@time plot_dist(files, varnames, varnames_print, components, Δt, nboxlength;
   interval, nplotstride)

println("Virtual satellite extraction done!")