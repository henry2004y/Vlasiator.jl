# ---
# title: Shared colorbar
# id: demo_2d_contour_share_colorbar
# date: 2023-10-16
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.3
# description: This demo shows how to share colorbar between subplots
# ---

This script shows how to share the same colorbar for multiple contour plots.
```julia
using Vlasiator, VlasiatorPyPlot

file = "bulk.vlsv"

meta = load(file)

pArgs = Vlasiator.set_args(meta, "fg_e", SI; normal=:none)

x1, x2 = Vlasiator.get_axis(pArgs)

norm = matplotlib.colors.CenteredNorm()

fig, axs = plt.subplots(1, 3, sharex=true, sharey=true, constrained_layout=true,
   figsize=(12,5))

data = Vlasiator.prep2d(meta, "Ehallx", 1)'
c = axs[1].pcolormesh(x1, x2, data; norm)

data = Vlasiator.prep2d(meta, "Ehally", 2)'
c = axs[2].pcolormesh(x1, x2, data; norm)

data = Vlasiator.prep2d(meta, "Ehallz", 3)'
c = axs[3].pcolormesh(x1, x2, data; norm)

for ax in axs
   ax.axis("scaled")
end

for i in 1:3
   axs[i].set_xlabel("x", fontsize=14)
end
for i in 1:1
   axs[i].set_ylabel("z", fontsize=14)
end

axs[1].set_title(L"$E_{hall,x}$", fontsize=16)
axs[2].set_title(L"$E_{hall,y}$", fontsize=16)
axs[3].set_title(L"$E_{hall,z}$", fontsize=16)

fig.colorbar(c, ax=axs[:])

savefig("Ehall.png", bbox_inches="tight")
```