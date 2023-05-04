# ---
# title: Vector components
# id: demo_vector_comps
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to plot vector components
# ---

This demo shows how to plot the components of a vector from one snapshot.
```julia
using Vlasiator, VlasiatorPyPlot, LaTeXStrings, Printf

function main()
   file = "/wrk/group/spacephysics/vlasiator/2D/ABA/bulk/bulk.0001000.vlsv"

   axisunit = EARTH
   colorscale = Linear
   cmap = matplotlib.cm.RdBu

   var = "B_vol"

   #####
   meta = load(file)

   pArgs = Vlasiator.set_args(meta, var, axisunit; normal=:none)

   x, y = Vlasiator.get_axis(pArgs)

   B = meta["B_vol"] # [T]
   bx = @views B[1,:].*1e9
   bx = reshape(bx, meta.ncells[1], meta.ncells[2])
   by = @views B[2,:].*1e9
   by = reshape(by, meta.ncells[1], meta.ncells[2])
   bz = @views B[3,:].*1e9
   bz = reshape(bz, meta.ncells[1], meta.ncells[2])

   V = (bx, by, bz)
   V_str = ("Bx", "By", "Bz")

   # Symmetric range for diverging colormap
   vmin = minimum(minimum.(V))
   vmax = -vmin

   cnorm1, cticks1 = set_colorbar(colorscale, vmin, vmax)

   fig, axs = subplots(1,3,
      figsize=(10, 4), sharex=true, sharey=true, constrained_layout=true)

   c1 = axs[1].pcolormesh(x, y, V[1]'; norm=cnorm1, cmap)

   for i in eachindex(axs)[2:end]
      axs[i].pcolormesh(x, y, V[i]'; norm=cnorm1, cmap)
   end

   for (i, ax) in enumerate(axs)
      ax.set_aspect("equal")
      ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
      ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
      ax.set_title(V_str[i])
      ax.set_xlabel(L"x [$R_E$]")
   end

   axs[1].set_ylabel(L"y [$R_E$]")

   # One colorbar for all subplots
   cb1 = colorbar(c1; ax=axs, ticks=cticks1, fraction=0.046, pad=0.04)
   cb1.ax.set_ylabel("[nT]")
   cb1.outline.set_linewidth(1.0)

   str_title = @sprintf "B at t = %4.1fs" meta.time

   fig.suptitle(str_title, fontsize="xx-large")

   savefig("Bvector.png", bbox_inches="tight")
end

main()
```