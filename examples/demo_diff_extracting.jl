# Sample script for extracting differences between FsGrid and DCCRG variables from 3D data.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot
import Vlasiator: RE

function main()
   directory = "./"
   file = "bulk1.0001284.vlsv"
   # Box range for variation diff extraction
   extent = [0.0, 12.0, -7.0, 7.0] # default [-Inf, Inf, -Inf, Inf]

   normal = :y
   origin = 0.0
   axisunit = EARTH
   comp = 2

   meta = load(directory*file)

   pArgs = Vlasiator.set_args(meta, "fg_e", axisunit; normal, origin)
   e_fg = Vlasiator.prep2dslice(meta, "fg_e", normal, comp, pArgs)

   pArgs = Vlasiator.set_args(meta, "vg_e_vol", axisunit; normal, origin)
   e_vg = Vlasiator.prep2dslice(meta, "vg_e_vol", normal, comp, pArgs)

   x1, x2 = Vlasiator.get_axis(pArgs)

   range1, range2 =
      if axisunit == EARTH
         searchsortedfirst(x1, extent[1]):searchsortedlast(x1, extent[2]),
         searchsortedfirst(x2, extent[3]):searchsortedlast(x2, extent[4])
      else
         searchsortedfirst(x1, extent[1]):searchsortedlast(x1, extent[2]),
         searchsortedfirst(x2, extent[3]):searchsortedlast(x2, extent[4])
      end

   ## Visualization
   fig, ax = plt.subplots()

   diff = (e_fg[range1,range2] - e_vg[range1,range2])'

   c = ax.pcolormesh(x1[range1], x2[range2], diff,
      norm=matplotlib.colors.CenteredNorm(),
      cmap=matplotlib.cm.RdBu_r)

   if comp == 1
      fig.colorbar(c; ax, label="E_fg_x - E_vg_x [V/m]")
   elseif comp == 2
      fig.colorbar(c; ax, label="E_fg_y - E_vg_y [V/m]")
   elseif comp == 3
      fig.colorbar(c; ax, label="E_fg_z - E_vg_z [V/m]")
   end

   ax.set_title(file)

   #savefig("test_diff.png")
end

main()