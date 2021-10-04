# Sample postprocessing script for wave-like structure spatial-temporal distribution
# searching.
#
# Usage:
#   julia -t 4 demo_wave_search.jl
# or
#   JULIA_NUM_THREADS=4 julia demo_wave_search.jl
#
# Procedures:
# 1. Extract variables in all cells from all snapshots.
# 2. For each cell and time interval, count time-series local peaks.
# 3. Plot the peak occurrence frequencies across the whole domain as an indicator for waves.
#
# Note:
# 1. When dealing with multiple variables, it is recommended to handle one variable at a
# time through the whole process due to memory considerations.
# 2. It assumes uniform sampling in time.
#
# Hongyang Zhou, hyzhou@umich.edu

using Glob, Vlasiator, PyPlot, LaTeXStrings

"Extract time series variable"
function extract_var(files, ncells, varname, component=0)
   nFiles = length(files)
   var = zeros(Float32, ncells[1], ncells[2], nFiles)

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
function countpeaks(y, n)
   minmaxs = Int[]
   for i in 2:length(y)-1
      if y[i+1] < y[i] > y[i-1] || y[i+1] > y[i] < y[i-1]
         push!(minmaxs, i)
      end
   end
   nCounts = zeros(Int, length(y)-n+1)
   nCounts[1] = count(i->(1 ≤ i ≤ n), minmaxs)
   for i in 1:length(y)-n
      nCounts[i+1] = nCounts[i]
      if i ∈ minmaxs
         nCounts[i+1] -= 1
      end
      if i+n-1 ∈ minmaxs
         nCounts[i+1] += 1
      end
   end
   nCounts
end

"Check wave-like occurrence frequencies within box length `n` of output interval `dt`."
function checkwaves_sma(var, dt=0.5, n::Int=size(var,3))
   nPeaks = zeros(Int, size(var,3)-n+1, size(var,1), size(var,2))

   Threads.@threads for j in axes(var, 2)
      for i in axes(var, 1)
         var_series = @view var[i,j,:]
         nPeaks[:,i,j] = countpeaks(var_series, n)
      end
   end
   nPeaks ./ (n*dt)
end

function plot_dist(files, varnames, varnames_print, components, Δt, nboxlength)
   @assert nboxlength ≥ 3 && isodd(nboxlength) "Expect odd box length ≥ 3!"
   if (local nFiles = length(files)) < nboxlength
      @warn "Set moving box length to the number of files..."
      nboxlength = nFiles
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

   for i in eachindex(varnames)
      # Obtain time series data
      var = extract_var(files, ncells, varnames[i], components[i])
      # Count local peak occuring frequencies at each location
      fPeaks = checkwaves_sma(var, Δt, nboxlength)
      for it in axes(fPeaks,1) # Iterate over time
         outname = "../out/spatial_perturbation_distribution_"*
            "$(lowercase(varnames_print[i]))_$(lpad(it, 3, '0')).png"
         isfile(outname) && continue
         ## Visualization
         im = ax.pcolormesh(y, x, fPeaks[it,:,:], shading="auto")
         ax.set_title("$(varnames_print[i]) Perturbation Detection, "*
            "t = $(round(tStart+(it-1)*Δt, digits=1)) ~ "*
            "$(round(tStart+(it+nboxlength-1)*Δt, digits=1))s";
            fontsize, fontweight="bold")
         ax.set_xlabel(L"y [$R_E$]"; fontsize, weight="black")
         ax.set_ylabel(L"x [$R_E$]"; fontsize, weight="black")
         ax.set_aspect("equal")
         ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
         ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
         ax.grid(true, color="grey", linestyle="-")

         cb = fig.colorbar(im; ax)
         cb.ax.set_ylabel("Frequency of local peak occurrence, [#/s]"; fontsize)

         savefig(outname, bbox_inches="tight")
         cla()
         cb.remove()
      end
   end
end

##### Main

varnames = ["proton/vg_rho", "vg_pressure", "proton/vg_v", "proton/vg_v", "vg_b_vol",
   "vg_e_vol", "vg_e_vol"]
varnames_print = ["Density", "Thermal Pressure", "Vx", "Vy", "Bz", "Ex", "Ey"]
components = [0, 0, 1, 2, 3, 1, 2] # 0: scalar; 1: x, 2: y, 3: z
Δt = 0.5                           # output time interval
nboxlength = 201                   # moving box average length
dir = "../run_rho2_bz-5_timevarying_startfrom300s" # data directory

files = glob("bulk*.vlsv", dir)

println("Total number of snapshots: $(length(files))")
println("Running with $(Threads.nthreads()) threads...")

@time plot_dist(files, varnames, varnames_print, components, Δt, nboxlength)

println("Virtual satellite extraction done!")