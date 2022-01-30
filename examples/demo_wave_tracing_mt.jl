# Script for tracing waves and plotting the dispersion relation.
#
# On a equatorial plane, B ∥ ẑ, we can choose an arbitrary line in-plane.
# On a meridional plane, B is in-plane, we can find a local line region ∥ B and ⟂ B.
#
# Currently only working on a equatorial plane.
# Usage:
#   julia -t 4 demo_wave_tracing_mt.jl
# or
#   JULIA_NUM_THREADS=4 julia demo_wave_tracing_mt.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator
using Vlasiator: qᵢ, μ₀, c, mᵢ, ϵ₀, RE
using Glob, DSP, FFTW, ImageFiltering, Interpolations
using Statistics: mean
using LinearAlgebra
using PyPlot

## Types

struct Variables
   varnames::Vector{String}
   varnames_print::Vector{String}
   components::Vector{Int}
end

## Methods

ispolar(meta::MetaVLSV) = findfirst(==(1), meta.ncells) == 2

"Extract `component` of variable `varname` at `cellids` from `files`."
function extract_var(files, varname, cellids, distances, component=0)
   sample_loc = range(distances[1], distances[end], length=length(distances))

   var = zeros(length(distances), length(files))

   Threads.@threads for i in eachindex(files)
      meta = load(files[i])
      if component == 0
         var_line = readvariable(meta, varname, cellids)[:]
      else
         var_line = readvariable(meta, varname, cellids)[component,:]
      end
      #TODO: do we need high order interpolations?
      interp_linear = LinearInterpolation(distances, var_line)
      var_line_resample = interp_linear.(sample_loc)
      var_line_smooth = imfilter(var_line_resample, Kernel.gaussian((3,)))

      var[:,i] = var_line_smooth
   end
   var
end

"CFL constrained normalized frequency."
dispersion_CFL(k, dx, dt, di, ωci) = dx/dt * abs(k) /(di * ωci)

"Normalized frequency of fast magnetosonic waves along angle `θ` with Doppler shift."
function dispersion_fast_perp(k, θ, vS, vA, v, di, ωci)
   ω = zeros(length(k))

   turnindex = findfirst(>=(0), k)

   vbulkpar = v[1]*cos(θ) + v[2]*sin(θ)

   dv1 =  √(vS^2 + vA^2) + vbulkpar # propagate along +θ direction
   dv2 = -√(vS^2 + vA^2) + vbulkpar # propagate along -θ direction

   if dv1 < 0; dv1 = 0.0; end
   if dv2 > 0; dv2 = 0.0; end

   for i in 1:turnindex-1
      ω[i] = dv2*k[i] /(di * ωci)
   end

   for i in turnindex:length(k)
      ω[i] = dv1*k[i] /(di * ωci)
   end
   ω
end

"Normalized frequency of bulk flow along tilted angle `θ`."
function dispersion_bulk_flow(k, θ, v, di, ωci)
   ω = zeros(length(k))
   turnindex = findfirst(>=(0), k)
   vbulkpar = v[1]*cos(θ) + v[2]*sin(θ)
   irange = vbulkpar > 0 ? (turnindex:length(k)) : (1:turnindex)
   for i in irange # otherwise 0
      ω[i] = vbulkpar*k[i] /(di * ωci)
   end
   ω
end

"Return the index in sorted `vec` with value closest to `x`."
function searchsortednearest(vec, x)
   idx = searchsortedfirst(vec, x)
   if idx == 1
      return idx
   elseif idx > length(vec)
      return length(vec)
   elseif vec[idx] == x
      return idx
   elseif abs(vec[idx]-x) < abs(vec[idx-1]-x)
      return idx
   else
      return idx-1
   end
end

"Obtain fast mode and bulk speed along line `angle` w.r.t. x-axis for cell ID `cid` in
`meta`."
function getCharacteristicSpeeds(meta, cid, angle)
   n = readvariable(meta, "proton/vg_rho", cid)
   B = readvariable(meta, "vg_b_vol", cid)
   p = readvariable(meta, "vg_pressure", cid)
   v = readvariable(meta, "proton/vg_v", cid)

   Bmag = norm.(eachcol(B))

   vA  = @. Bmag / √(μ₀ * n * mᵢ)        # Alfven speed, [m/s]
   vS  = @. √(γ * p / (n * mᵢ))          # sonic speed, [m/s]

   vFast = @. √(vA^2 + vS^2)
   vBulk = @. v[1,:]*cos(angle) + v[2,:]*sin(angle)

   vFast, vBulk
end

"Trace along line with `angle` w.r.t. x-axis at possible wave speeds."
function tracewave(xstart, angle, files, dtfile, t, cellids, distances)

   x1 = zeros(length(t)); x1[1] = xstart
   x2 = copy(x1)
   x3 = copy(x1)

   tfile = let
      tfilefirst = getproperty(load(files[1]), :time)
      tfilelast = getproperty(load(files[end]), :time)
      tfilefirst:dtfile:tfilelast
   end

   ifile = 1 # file index tracker

   dt = t[2] - t[1] # timestep

   for it in eachindex(t)[1:end-1]
      # Find the closest output saving time to t[it], 0th order interpolation
      for i = ifile:length(tfile)
         if abs(tfile[i] - t[it]) < 0.5*dtfile
            ifile = i
            break
         end
      end
      meta = load(files[ifile])
      # Find the closest cell to the wave location
      cid = let
         c1 = searchsortednearest(distances, x1[it])
         c2 = searchsortednearest(distances, x2[it])
         c3 = searchsortednearest(distances, x3[it])
         cellids[c1], cellids[c2], cellids[c3]
      end
      vFast, vBulk = getCharacteristicSpeeds(meta, cid, angle)
      x1[it+1] = x1[it] + dt * vBulk[1]
      x2[it+1] = x2[it] + dt * (vBulk[2] + vFast[2])
      x3[it+1] = x3[it] + dt * (vBulk[3] - vFast[3])
   end

   x1, x2, x3
end

"Evaluate the average quantities at `cellids` at the middle of `files`."
function estimate_meanstates(files, cellids)
   # Select the snapshot in the middle
   nfile = length(files)
   meta = load(files[nfile÷2+1])

   n = readvariable(meta, "proton/vg_rho", cellids)
   v = readvariable(meta, "proton/vg_v", cellids)
   p = readvariable(meta, "vg_pressure", cellids)

   if hasvariable(meta, "vg_b_vol")
      B = readvariable(meta, "vg_b_vol", cellids)
   elseif hasvariable(meta, "fg_b")
      B = readvariable(meta, "fg_b", cellids)
   else
      B = readvariable(meta, "b", cellids)
   end

   vperp = @view v[1:2,:]
   vpar = @view v[3,:]

   # Obtain average states
   n̄ = mean(n)
   p̄ = mean(p)
   v̄par = mean(vpar)
   v̄perp = @views [mean(vperp[1,:]), mean(vperp[2,:])]

   # Characteristic parameters
   Bnorm = @views abs(mean(B[3,:]))
   di  = √(mᵢ*ϵ₀/(n̄))*c/qᵢ               # ion inertial length, [m]
   ωci = qᵢ*Bnorm/mᵢ                     # [/s]
   v̄A  = Bnorm / √(μ₀ * n̄ * mᵢ)          # Alfven speed, [m/s]
   v̄S  = √(γ * p̄ / (n̄ * mᵢ))             # sonic speed, [m/s]

   println("--------------------------------------------------")
   println("* Average states along the line at the middle snapshot")
   println("Density               : ", rpad(round(n̄/1e6; digits=2), 8), "amu/cc")
   println("Pressure              : ", rpad(round(p̄*1e9; digits=3), 8), "nPa")
   println("Parallel velocity     : ", rpad(round(v̄par/1e3; digits=2), 8), "km/s")
   println("Perpendicular velocity: ", rpad(round.(v̄perp/1e3; digits=2), 8), "km/s")
   println("Flow angle            : ", rpad(round(atand(v̄perp[2], v̄perp[1]); digits=2), 8),
      "degrees")
   println("Ion inertial length   : ", rpad(round(di/1e3; digits=2), 8), "km")
   println("Gyrofrequency         : ", rpad(round(ωci; digits=2), 8), "Hz")
   println("Alfven speed          : ", rpad(round(v̄A/1e3; digits=2), 8), "km/s")
   println("Sonic speed           : ", rpad(round(v̄S/1e3; digits=2), 8), "km/s")
   println("--------------------------------------------------")

   di, ωci, v̄A, v̄S, v̄perp
end

"Plot the process of wave checks."
function plot_dispersion(files, vars, cellids, distances, coords, meanstates, dtfile, Δt)

   # Parameters
   nfile = length(files)
   npoints = length(cellids)
   nt = nfile ÷ 2 + 1

   varnames = vars.varnames
   varnames_print = vars.varnames_print
   components = vars.components

   di, ωci, v̄A, v̄S, v̄perp = meanstates

   tfile1st = load(files[1]).time
   # Output timestamps
   t = [dtfile * ifile + tfile1st for ifile in 0:nfile-1]
   # Selected line tilted angle ∈ [-π, π]
   θ = atan(coords[2,end] - coords[2,1], coords[1,end] - coords[1,1])
   # Sample width, [m]
   dx = norm(coords[:,end] .- coords[:,1]) /(npoints - 1)
   println("spatial resolution: ", round(dx/1e3; digits=2), " km")

   # Trace wave along the line in a space-time domain
   twave = let
      tstart = 400.0
      tend = 430.0
      tstart:200Δt:tend
   end
   xwave = tracewave(distances[npoints÷2+1], θ, files, dtfile, twave, cellids, distances)

   twave2 = let
      tstart = 600.0
      tend = 630.0
      tstart:200Δt:tend
   end
   xwave2 = tracewave(distances[npoints÷2+1], θ, files, dtfile, twave2, cellids, distances)

   # Dispersion plotting ranges
   kmin = -π / dx * di         # minimum wave number
   kmax =  π / dx * di         # maximum wave number
   ωmin = 0                    # minimum angular frequency
   ωmax = π / dtfile / ωci     # maximum angular frequency

   # Only the 1st quadrature
   krange = range(kmin, kmax, length=npoints)
   ωrange = range(ωmin, ωmax, length=nt)

   axisunit = RE

   # Precalculated lines
   ωCFL = dispersion_CFL.(krange, dx, Δt, di, ωci)
   ωfast = dispersion_fast_perp(krange, θ, v̄S, v̄A, v̄perp, di, ωci)
   ωbulk = dispersion_bulk_flow(krange, θ, v̄perp, di, ωci)
   # Window filtering for avoiding spectral leakage
   window = hanning(npoints) * hanning(nfile)'

   meta = load(files[end])

   for i in eachindex(varnames)
      println("variable name: ", varnames[i])
      var = extract_var(files, varnames[i], cellids, distances, components[i])

      # 2DFFT
      F̃ = window .* var |> fft |> fftshift

      # Visualization
      fig = figure(figsize=(12,12), constrained_layout=false)
      ax = [subplot(221), subplot(223), subplot(222), subplot(224, projection="3d")]

      dispersion = reverse!(abs.(F̃.*F̃)[:, end-nt+1:end]', dims=1)
      im1 = ax[1].pcolormesh(krange, ωrange, dispersion, norm=matplotlib.colors.LogNorm())

      ax[1].plot([krange[1], 0.0, krange[end]], [ωCFL[1], 0.0, ωCFL[end]], "--",
         linewidth=1.0, color="k", label="CFL Condition")
      ax[1].plot(krange, ωfast, "--",
         linewidth=1.2, color="#d62728", label="Fast Mode")
      ax[1].plot(krange, ωbulk, "--",
         linewidth=1.2, color="#9467bd", label="Flow Speed")

      ax[1].set_xlim(kmin, kmax)
      ax[1].set_ylim(ωmin, ωmax)

      cb = colorbar(im1; ax=ax[1])
      cb.ax.tick_params(direction="in")
      ax[1].legend(;fontsize="x-small")
      ax[1].set_xlabel(L"$k_\perp \cdot d_i$")
      ax[1].set_ylabel(L"$\omega/\Omega_{ci}$")
      ax[1].set_title(L"$k_\perp$ angle w.r.t. x = %$θ")

      im2 = ax[2].pcolormesh((distances .+ coords[1,1])./RE, t, var')

      ax[2].plot((xwave[1] .+ coords[1,1])./RE, twave, ".--",
         color="#d62728", label=L"$V_{bulk}$")
      ax[2].plot((xwave[2] .+ coords[1,1])./RE, twave, ".--",
         color="#9467bd",  label=L"$V_{bulk} + V_{fast}$")
      ax[2].plot((xwave[3] .+ coords[1,1])./RE, twave, ".--",
         color="#ff7f0e", label=L"$V_{bulk} - V_{fast}$")

      ax[2].plot((xwave2[1] .+ coords[1,1])./RE, twave2, ".--",
         color="#d62728")
      ax[2].plot((xwave2[2] .+ coords[1,1])./RE, twave2, ".--",
         color="#9467bd")
      ax[2].plot((xwave2[3] .+ coords[1,1])./RE, twave2, ".--",
         color="#ff7f0e")

      cb = colorbar(im2; ax=ax[2])
      cb.ax.tick_params(direction="in")

      ax[2].set_xlim(coords[1,1]/RE, coords[1,end]/RE)

      ax[2].legend(loc="upper center", bbox_to_anchor=(0.5, -0.13),
         fancybox=true, shadow=true, ncol=3)
      ax[2].set_xlabel(L"x [$R_E$]")
      ax[2].set_ylabel(L"time [s]")
      ax[2].set_title("$(varnames_print[i])_$(components[i])")

      pArgs = Vlasiator.set_args(meta, varnames[i], axisunit; normal=:none)
      x, y = Vlasiator.get_axis(pArgs)
      data = Vlasiator.prep2d(meta, varnames[i], components[i])'
      cnorm, cticks = Vlasiator.set_colorbar(Vlasiator.Linear, -Inf, Inf, data)
      cmesh = ax[3].pcolormesh(x, y, data, norm=cnorm)

      ax[3].set_xlabel(pArgs.strx)
      ax[3].set_ylabel(pArgs.stry)
      ax[3].set_aspect(1)

      cb = colorbar(cmesh; ax=ax[3], ticks=cticks, fraction=0.046, pad=0.04)
      cb.ax.set_ylabel(pArgs.cb_title)
      cb.ax.tick_params(direction="in")

      @views ax[3].scatter(coords[1,:]./RE, coords[2,:]./RE; s=0.2, color="k")

      xCoord = (distances .+ coords[1,1])./RE
      # meshgrid
      X = [x for _ in t, x in xCoord]
      Y = [y for y in t, _ in xCoord]

      ax[4].view_init(elev=40., azim=-30.)
      ax[4].plot_surface(X, Y, var', cmap=matplotlib.cm.turbo, antialiased=false)

      ax[4].set_xlabel(L"x [$R_E$]"; fontsize="small")
      ax[4].set_ylabel("time [s]"; fontsize="small")
      #ax[4].set_zlabel(pArgs.cb_title)
      ax[4].tick_params(labelsize="small")

      outname = "dispersion_$(varnames_print[i])_$(components[i]).png"
      savefig(joinpath(outdir, outname), bbox_inches="tight")
      close(fig)
   end

end

## Main

const outdir = "../out/"
const γ = 5 / 3

xStart = [10.0, 0.0, 0.0].*RE
xEnd   = [13.3, 0.0, 0.0].*RE

varnames = ["proton/vg_rho", "vg_b_vol", "vg_e_vol", "vg_pressure"]
varnames_print = ["rho", "b", "e", "p"]
components = [0, 3, 1, 0] # 0 for scalar, 1-3 for vector components

dir = "../run_rho2_bz-5_timevarying_startfrom300s"
files = glob("bulk.*.vlsv", dir)

vars = Variables(varnames, varnames_print, components)

dtfile = 0.5                          # output interval, [s]
Δt   = 0.0147176                      # discrete time step from runlog, [s]

meta = load(files[1])

if ispolar(meta) # polar plane
   @error "not implemented!"
end

cellids, distances, coords = getcellinline(meta, xStart, xEnd)

meanstates = estimate_meanstates(files, cellids)

println("number of extracted points: ", length(cellids))
println("xStart: ", xStart)
println("xEnd: ", xEnd)
tbegin = load(files[1]).time
tend = load(files[end]).time
println("time from $tbegin to $tend s")

@time plot_dispersion(files, vars, cellids, distances, coords, meanstates, dtfile, Δt)