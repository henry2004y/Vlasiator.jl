# ---
# title: Poynting flux
# id: demo_poynting
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to extract poynting flux on a subdomain
# ---

This demo shows how to extract the Poynting flux on a subdomain from 2D and plot.

Usage:
```shell
julia -t 4 demo_poynting_mt_pyplot.jl
```
or
```shell
JULIA_NUM_THREADS=4 julia demo_poynting_mt_pyplot.jl
```

```julia
using Statistics: mean, normalize
using LinearAlgebra: ×, ⋅
using Vlasiator
using Vlasiator: μ₀, prep2d, prep2dslice
using DSP, Polyester, Printf, PyPlot

"""
    extract_EM(files, range1, range2, pArgs, ndim=2; normal=:y, verbose=true)

Extract time-series of EM fields from `files` within `range1` in the 1st dim and `range1` in
the 2nd dim.
"""
function extract_EM(files, range1, range2, pArgs, ndim=2; normal=:y, verbose=true)
   @assert ndim ∈ (2,3) "Only support extracting EM fields from 2D/3D runs."
   nfile = length(files)
   e = zeros(Float32, 3, length(range1), length(range2), nfile)
   b = zeros(Float32, 3, length(range1), length(range2), nfile)

   if ndim == 2
      @batch for i in eachindex(files)
         verbose && println("i = $i/$nfile, file = $(files[i])")
         meta = load(files[i])

         e[1,:,:,i] = prep2d(meta, "vg_e_vol", 1)[range1, range2]
         e[2,:,:,i] = prep2d(meta, "vg_e_vol", 2)[range1, range2]
         e[3,:,:,i] = prep2d(meta, "vg_e_vol", 3)[range1, range2]
         b[1,:,:,i] = prep2d(meta, "vg_b_vol", 1)[range1, range2]
         b[2,:,:,i] = prep2d(meta, "vg_b_vol", 2)[range1, range2]
         b[3,:,:,i] = prep2d(meta, "vg_b_vol", 3)[range1, range2]
      end
   else
      @batch for i in eachindex(files)
         verbose && println("i = $i/$nfile, file = $(files[i])")
         meta = load(files[i])

         e[1,:,:,i] = prep2dslice(meta, "vg_e_vol", normal, 1, pArgs)[range1, range2]
         e[2,:,:,i] = prep2dslice(meta, "vg_e_vol", normal, 2, pArgs)[range1, range2]
         e[3,:,:,i] = prep2dslice(meta, "vg_e_vol", normal, 3, pArgs)[range1, range2]
         b[1,:,:,i] = prep2dslice(meta, "vg_b_vol", normal, 1, pArgs)[range1, range2]
         b[2,:,:,i] = prep2dslice(meta, "vg_b_vol", normal, 2, pArgs)[range1, range2]
         b[3,:,:,i] = prep2dslice(meta, "vg_b_vol", normal, 3, pArgs)[range1, range2]
      end
   end
   e, b
end

"Moving box average of vector `g` with box size `n`."
function moving_average(g::AbstractVector{<:AbstractFloat}, n::Int)
   nbackward = div(n,2)
   nforward = isodd(n) ? div(n,2) : div(n,2) - 1
   len = length(g)
   g_avg = similar(g)
   @inbounds @batch for i = 1:len
       lo = max(1, i - nbackward)
       hi = min(len, i + nforward)
       g_avg[i] = mean(@view g[lo:hi])
   end
   g_avg
end

"Extract the perturbation from vector `v` with moving box size `nbox`."
function detrend(v::AbstractVector{<:AbstractFloat}; nbox=length(v)÷12)
   v̄ = moving_average(v, nbox)
   dv = v .- v̄
end

"Extract the perturbation from 4D array `v` along the last dim with moving box size `nbox`."
function detrend(v::AbstractArray{<:AbstractFloat}; nbox=size(v,4)÷12)
   v̄ = similar(v)
   for idim = 1:3
      for j in axes(v, 3), i in axes(v, 2)
         v̄[idim,i,j,:] = @views moving_average(v[idim,i,j,:], nbox)
      end
   end
   dv = v .- v̄
   dv, v̄
end

"Band pass filter for vector field `var`."
function band_pass_filter(var, responsetype, designmethod)

   filter = digitalfilter(responsetype, designmethod)

   var_filtered = zeros(eltype(var), size(var,4), 3, size(var,2), size(var,3))

   for idim = 1:3
      for j in axes(var, 3), i in axes(var, 2)
         var_filtered[:,idim,i,j] = filtfilt(filter, var[idim,i,j,:])
      end
   end
   var_filtered
end

"""
    calc_poynting(de, db, b̄, it)

Return full Poynting vector, also parallel and perpendicular Poynting vector components to
the magnetic field. The perpendicular components are in-plane.
"""
function calc_poynting(de, db, b̄, it)

   s = zeros(3, size(de,3), size(de,4))
   # parallel and perpendicular Poynting vector magnitudes
   s_par  = zeros(size(de,3), size(de,4))
   s_perp = zeros(size(de,3), size(de,4))

   @views for j in axes(de,4), i in axes(de,3)
      # Poynting vector = dE × dB
      s[:,i,j] = de[it,:,i,j] × db[it,:,i,j] / μ₀

      # Transform into parallel and perpendicular direction w.r.t. B
      b̂ = normalize(b̄[:,i,j,it])
      s_par[i,j] = s[:,i,j] ⋅ b̂
      #s_perp[i,j] = norm(s[:,i,j] .- s_par[i,j] .* b̂)
      # in-plane perpendicular component only
      s_perp[i,j] = sqrt((s[1,i,j] - s_par[i,j]*b̂[1])^2 + (s[3,i,j] - s_par[i,j]*b̂[3])^2)
   end

   s, s_par, s_perp
end

"Output figures of Poynting vectors for each snapshot."
function plot_poynting(de_filtered, db_filtered, b̄, x1, x2, frequency_range)

   nt = size(de_filtered, 1)

   norm1 = let
      sparmax = frequency_range == "high" ? 5e-7 : 1e-5
      sparmin = -sparmax
      levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(sparmin, sparmax)
      matplotlib.colors.BoundaryNorm(levels, ncolors=256, clip=true)
   end

   norm2 = let
      sperpmax = frequency_range == "high" ? 5e-7 : 1e-5
      sperpmin = 0.0
      levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(sperpmin, sperpmax)
      matplotlib.colors.BoundaryNorm(levels, ncolors=256, clip=true)
   end

   stride = 10 # number of strides for quivers

   s, s_par, s_perp = calc_poynting(de_filtered, db_filtered, b̄, 1)


   fig, axs = plt.subplots(1, 2; figsize=(11,6), sharex=true, sharey=true,
      constrained_layout=true)

   axs[1].set_title(L"S_\parallel"; fontsize="large")
   axs[2].set_title(L"S_\perp"; fontsize="large")

   for ax in axs
      ax.set_aspect("equal")
      ax.set_xlabel(L"x1 [R_E]"; fontsize="large")
      ax.set_ylabel(L"x2 [R_E]"; fontsize="large")
   end

   c1 = axs[1].pcolormesh(x1, x2, s_par';
      cmap=matplotlib.cm.RdBu_r,
      shading="nearest",
      norm=norm1,
      )

   cb1 = colorbar(c1; ax=axs[1], extend="both")
   cb1.ax.set_ylabel(L"[W/m^2]"; fontsize="large")

   c2 = axs[2].pcolormesh(x1, x2, s_perp';
      cmap=matplotlib.cm.turbo,
      shading="nearest",
      norm=norm2,
      )

   cb2 = colorbar(c2; ax=axs[2], extend="max")
   cb2.ax.set_ylabel(L"[W/m^2]"; fontsize="large")

   s1 = @views (s[1,:,:] ./ hypot.(s[1,:,:], s[3,:,:]))[1:stride:end, 1:stride:end]'
   s2 = @views (s[3,:,:] ./ hypot.(s[1,:,:], s[3,:,:]))[1:stride:end, 1:stride:end]'

   q = axs[1].quiver(x1[1:stride:end], x2[1:stride:end], s1, s2; color="k")

   b1 = @views (b̄[1,:,:,1] ./ hypot.(b̄[1,:,:,1], b̄[3,:,:,1]))[1:stride:end, 1:stride:end]'
   b2 = @views (b̄[3,:,:,1] ./ hypot.(b̄[1,:,:,1], b̄[3,:,:,1]))[1:stride:end, 1:stride:end]'

   qb = axs[1].quiver(x1[1:stride:end], x2[1:stride:end], b1, b2; color="tab:purple")

   if frequency_range == "high"
      fig.suptitle("t = $(t[1]) s, [0.1, 1.0] Hz";
         fontsize="xx-large")
   else
      fig.suptitle("t = $(t[1]) s, [$(responsetype.w1), $(responsetype.w2)] Hz";
         fontsize="xx-large")
   end

   for it = 1:nt
      @info "it = $it"
      local s, s_par, s_perp
      outname = "poynting/poynting_$(lpad(it, 4, '0'))_band$frequency_range.png"
      isfile(outname) && continue
      s, s_par, s_perp = calc_poynting(de_filtered, db_filtered, b̄, it)

      c1.set_array(s_par')
      c2.set_array(s_perp')

      ŝx = @views @. s[1,:,:] / hypot(s[1,:,:], s[3,:,:])
      ŝz = @views @. s[3,:,:] / hypot(s[1,:,:], s[3,:,:])

      q.set_UVC(ŝx[1:stride:end, 1:stride:end]', ŝz[1:stride:end, 1:stride:end]')

      b̂x = @views @. b̄[1,:,:,it] / hypot(b̄[1,:,:,it], b̄[3,:,:,it])
      b̂z = @views @. b̄[3,:,:,it] / hypot(b̄[1,:,:,it], b̄[3,:,:,it])

      qb.set_UVC(b̂x[1:stride:end, 1:stride:end]', b̂z[1:stride:end, 1:stride:end]')

      if frequency_range == "high"
         fig.suptitle("t = $(t[it]) s, [0.1, 1.0] Hz";
            fontsize="xx-large")
      else
         fig.suptitle("t = $(t[it]) s, [$(responsetype.w1), $(responsetype.w2)] Hz";
            fontsize="xx-large")
      end
      savefig(outname; dpi=300, bbox_inches="tight")
   end
end

########## Main

files = filter(contains(r"^bulk.*\.vlsv$"), readdir())
nfile = length(files)

const frequency_range = "low" # filtered frequency range ∈ ("low", "high")
const extent = [6., 16., -7., 7.] # [RE], default: [-Inf32, Inf32, -Inf32, Inf32]
const normal = :none # plane normal direction for 3D data, (:none, :x, :y, :z)
const fs = 2.0 # sampling frequency, [Hz]

# WARNING: t may not be exact due to round-off errors in output time stamps!
ndim, pArgs, t, range1, range2 = let meta = load(files[1])
   pArgs = Vlasiator.set_args(meta, "vg_e_vol", EARTH; normal=:none)
   x1, x2 = Vlasiator.get_axis(pArgs)

   tstart = meta.time
   tend = load(files[end]).time

   ndims(meta), pArgs,
   range(round(tstart,digits=1), round(tend, digits=1), length=nfile),
   searchsortedfirst(x1, extent[1]):searchsortedlast(x1, extent[2]),
   searchsortedfirst(x2, extent[3]):searchsortedlast(x2, extent[4])
end

@assert ndim == 3 "3D not working yet!"

e, b = extract_EM(files, range1, range2, pArgs, ndim; normal, verbose=true)

designmethod = Butterworth(5)
if frequency_range == "high"
   responsetype = Highpass(0.1; fs)
   nbox = 40
elseif frequency_range == "low"
   responsetype = Bandpass(0.02, 0.067; fs)
   nbox = 200
end

de, _ = detrend(e; nbox)
db, b̄ = detrend(b; nbox)

de_filtered = band_pass_filter(de, responsetype, designmethod)
db_filtered = band_pass_filter(db, responsetype, designmethod)

x1 = range(extent[1], extent[2], length=length(range1))
x2 = range(extent[3], extent[4], length=length(range2))

plot_poynting(de_filtered, db_filtered, b̄, x1, x2, frequency_range)
```