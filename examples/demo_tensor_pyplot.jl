# Sample postprocessing script for plotting a tensor from one snapshot.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot, LaTeXStrings, Printf

file = "/home/hongyang/Vlasiator/test/equatorial/timevarying/bulk.0001582.vlsv"

axisunit = RE
colorscale = Linear
cmap = matplotlib.cm.turbo

vars = ["proton/vg_ptensor_diagonal", "proton/vg_ptensor_offdiagonal"]

#################

meta = load(file)

pArgs1 = Vlasiator.set_args(meta, vars[1], axisunit; normal=:none)
pArgs2 = Vlasiator.set_args(meta, vars[2], axisunit; normal=:none)

x, y = Vlasiator.get_axis(pArgs1)
# [nPa]
pxx = Vlasiator.prep2d(meta, vars[1], :x)' .* 1e9 
pyy = Vlasiator.prep2d(meta, vars[1], :y)' .* 1e9
pzz = Vlasiator.prep2d(meta, vars[1], :z)' .* 1e9
pyz = Vlasiator.prep2d(meta, vars[2], :1)' .* 1e9
pxz = Vlasiator.prep2d(meta, vars[2], :2)' .* 1e9
pxy = Vlasiator.prep2d(meta, vars[2], :3)' .* 1e9

P = (pxx, pyy, pzz, pyz, pxz, pxy)
P_str = ("Pxx", "Pyy", "Pzz", "Pyz", "Pxz", "Pxy")

vmin = minimum(minimum.(P))
vmax = maximum(maximum.(P))

cnorm1, cticks1 = Vlasiator.set_colorbar(colorscale, vmin, vmax)

fig, axs = subplots(3,2,
   figsize=(6, 12), sharex=true, sharey=true, constrained_layout=true)

c1 = axs[1].pcolormesh(x, y, pxx; norm=cnorm1, cmap, shading="nearest")

for i in eachindex(axs)[2:end]
   axs[i].pcolormesh(x, y, P[i]; norm=cnorm1, cmap, shading="nearest")
end

for (i, ax) in enumerate(axs)
   ax.set_aspect("equal")
   ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
   ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
   ax.set_title(P_str[i])
end

for ax in axs[3,:]
   ax.set_xlabel(L"x [$R_E$]")
end

for ax in axs[:,1]
   ax.set_ylabel(L"y [$R_E$]")
end

# One colorbar for all subplots
cb1 = colorbar(c1; ax=axs, ticks=cticks1, fraction=0.046, pad=0.04)
cb1.ax.set_ylabel("[nPa]")
cb1.outline.set_linewidth(1.0)

str_title = @sprintf "Pressure Tensor at t = %4.1fs" meta.time

fig.suptitle(str_title, fontsize="xx-large")

savefig("ptensor.png", bbox_inches="tight")