# Identification of EMIC and mirror modes in two frequency bands from single satellite data.
#
# EMIC: δB⟂ vs. δv_⟂
# Mirror: δB∥ vs. δn 
#
# The moving box window size is decided empirically.
#
# Hongyang Zhou, hyzhou@umich.edu

using Statistics: mean
using LinearAlgebra: ⋅, normalize!, norm
using DSP, PyPlot, Glob
using Vlasiator
using Vlasiator: μ₀, mᵢ, RE

"Extract variables at `loc` from VLSV `files`."
function extract_vars(files, loc)
   nfiles = length(files)

   t = zeros(Float32, nfiles)
   n = zeros(Float32, nfiles)
   v = zeros(Float32, 3, nfiles)
   b = zeros(Float32, 3, nfiles)
   e = zeros(Float32, 3, nfiles)

   id = load(files[1]) do meta
      getcell(meta, loc)
   end

   # Extract data from each frame
   for i = eachindex(files)
      meta = load(files[i])
      t[i] = meta.time
      n[i] = readvariable(meta, "proton/vg_rho", id)[1]
      v[:,i] = readvariable(meta, "proton/vg_v", id)
      b[:,i] = readvariable(meta, "vg_b_vol", id)
      e[:,i] = readvariable(meta, "vg_e_vol", id)
   end

   t, n, v, e, b
end

function moving_average(g::AbstractVector{<:AbstractFloat}, n::Int)
   nbackward = div(n,2) 
   nforward = isodd(n) ? nbackward : nbackward - 1
   len = length(g)
   g_avg = similar(g)
   @inbounds for i = 1:len
      lo = max(1, i - nbackward)
      hi = min(len, i + nforward)
      g_avg[i] = mean(@view g[lo:hi])
   end
   g_avg
end

function moving_average(g::AbstractMatrix{<:AbstractFloat}, n::Int)
   nbackward = div(n,2) 
   nforward = isodd(n) ? nbackward : nbackward - 1
   len = size(g,2)
   g_avg = similar(g)
   for icomp = axes(g,1)
      @inbounds for i = 1:len
         lo = max(1, i - nbackward)
         hi = min(len, i + nforward)
         g_avg[icomp,i] = mean(@view g[icomp, lo:hi])
      end
   end
   g_avg
end

function detrend(v::AbstractVector{<:AbstractFloat}; nbox=length(v)÷12)
   v̄ = moving_average(v, nbox)
   dv = v .- v̄
end

function detrend(v::AbstractMatrix{<:AbstractFloat}; nbox=size(v,2)÷12)
   v̄ = similar(v)
   for i in axes(v,1)
      v̄[i,:] = @views moving_average(v[i,:], nbox)
   end
   dv = v .- v̄
end

function align_yaxis(ax1, ax2)
   y_lims = [ax.get_ylim() for ax in [ax1, ax2]]

   # normalize both axes
   y_mags = ntuple(i -> y_lims[i][2] - y_lims[i][1], Val(2))
   y_lims_normalized = ntuple(i -> y_lims[i] ./ y_mags[i], Val(2))

   # find combined range
   y_new_lims_normalized = (min(y_lims_normalized[1][1],y_lims_normalized[2][1]),
      max(y_lims_normalized[1][2], y_lims_normalized[2][2]))

   # denormalize combined range to get new axes
   new_lim1 = y_new_lims_normalized .* y_mags[1]
   new_lim2 = y_new_lims_normalized .* y_mags[2]
   ax1.set_ylim(new_lim1)
   ax2.set_ylim(new_lim2)
   return
end

################# Main

files = glob("bulk*.vlsv", ".")

# virtual satellite location
loc = [12Vlasiator.RE, 0, 0]

nbox = 40 # moving box size, [# of indices]
fs = 2.0 # sampling rate, [Hz]

println("Number of files: $(length(files))")
println("Extracting location: $loc [m]")
println("Sampling rate: $fs [Hz]")
println("Moving box window size: $(nbox / fs) [s]")

t, n, v, e, b = extract_vars(files, loc)

bmag = [sqrt(b[1,i]^2 + b[2,i]^2 + b[3,i]^2) for i in eachindex(n)]

# Detrend before filtering to remove the lowest frequency changes
dn = detrend(n; nbox)
dbmag = detrend(bmag; nbox)

n̄ = moving_average(n, nbox)
b̄mag = moving_average(bmag, nbox)

v̄a = @. b̄mag / sqrt(μ₀ * mᵢ * n̄)

dv = detrend(v; nbox) # [m/s]
db = detrend(b; nbox) # [T]

# Decompose into parallel and perpendicular vector components
b̂₀ = moving_average(b, nbox)
for i in axes(b̂₀, 2)
   normalize!(@view b̂₀[:,i])
end

# δB∥
db_par = [db[:,i] ⋅ b̂₀[:,i] for i in axes(db, 2)]
# δB⟂
db_perp = [norm(db[:,i] .- db_par[i] .* b̂₀[:,i]) for i in axes(db, 2)]
# δv∥
dv_par = [dv[:,i] ⋅ b̂₀[:,i] for i in axes(dv, 2)]
# δv⟂
dv_perp = [norm(dv[:,i] .- dv_par[i] .* b̂₀[:,i]) for i in axes(dv, 2)]

designmethod = Butterworth(5)

responsetype = Highpass(0.1; fs)
dn_high = filtfilt(digitalfilter(responsetype, designmethod), dn)
dbpar_high = filtfilt(digitalfilter(responsetype, designmethod), db_par)
dbperp_high = filtfilt(digitalfilter(responsetype, designmethod), db_perp)
dvperp_high = filtfilt(digitalfilter(responsetype, designmethod), dv_perp)

responsetype = Bandpass(0.02, 0.067; fs)
dn_low = filtfilt(digitalfilter(responsetype, designmethod), dn)
dbpar_low = filtfilt(digitalfilter(responsetype, designmethod), db_par)

color1 = "tab:blue"
color2 = "tab:red"

fig, axs = plt.subplots(3, 1; figsize=(13, 9), sharex=true, constrained_layout=true)

fig.suptitle("Virtual satellite at $(string(round.(loc./RE, digits=2))) "*L"R_E";
   fontsize="x-large")

axs[1].plot(t, dn_low ./ n̄, color1, label=L"\delta n, 0.02-0.067\, Hz")
ax12 = axs[1].twinx()
ax12.plot(t, dbpar_low ./ b̄mag, color=color2, label=L"\delta B, 0.02-0.067\, Hz")

axs[2].plot(t, (dn_high ./ n̄), color1, label=L"\delta n, 0.1-1.0\, Hz")
ax22 = axs[2].twinx()
ax22.plot(t, (dbpar_high ./ b̄mag), color=color2, label=L"\delta B, 0.1-1.0\, Hz")

axs[3].plot(t, dvperp_high ./ v̄a, color1, label=L"\delta v_\perp, 0.1-1.0\, Hz")
ax32 = axs[3].twinx()
ax32.plot(t, dbperp_high ./ b̄mag, color=color2, label=L"\delta B, 0.1-1.0\, Hz")

axs[3].set_xlabel("time [s]"; fontsize="large")

n_str = (L"\delta n /n", L"\delta n / n", L"\delta v_\perp / V_A")

for (i, a) in enumerate(axs)
   a.hlines(0.0, t[1], t[end]; colors="k", linestyle="dashed", alpha=0.3)
   a.tick_params(axis="y", labelcolor=color1)
   a.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
   a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
   a.grid(true)
   a.set_ylabel(n_str[i], color=color1, fontsize=14)
end

axtwin = (ax12, ax22, ax32)
b_str = ( L"\delta B_\parallel / B_0", L"\delta B_\parallel / B_0", L"\delta B_\perp / B_0")

for (i, a) in enumerate(axtwin)
   a.tick_params(axis="y", labelcolor=color2)
   a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
   a.set_ylabel(b_str[i], color=color2, fontsize="large")
   align_yaxis(axs[i], a)
end

AnchoredText = matplotlib.offsetbox.AnchoredText

at = AnchoredText(
   "[0.02, 0.067] Hz band", prop=Dict("size"=>"medium"), frameon=true, loc="lower left")
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
axs[1].add_artist(at)

at = AnchoredText(
   "[0.1, 1.0] Hz band", prop=Dict("size"=>"medium"), frameon=true, loc="lower left")
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
axs[2].add_artist(at)

at = AnchoredText(
   "[0.1, 1.0] Hz band", prop=Dict("size"=>"medium"), frameon=true, loc="lower left")
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
axs[3].add_artist(at)

savefig("virtual_satellite_wave.png"; dpi=300)