# ---
# title: Field lines with customized seeds
# id: demo_2d_fieldlines
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to plot field lines with handpicked seeds.
# ---

This demo shows how to plot field lines with handpicked seeds.
```julia
using VlasiatorPyPlot, FieldTracer
using Vlasiator: RE # Earth radius, [m]

function main()
   file = "../../bulk.0000999.vlsv"
   meta = load(file)
   (;coordmin, coordmax, ncells) = meta
   pcolormesh(meta, "proton/vg_rho", colorscale=Linear)

   dim_ =
      if ncells[2] == 1
         (1,3)
      elseif ncells[3] == 1
         (1,2)
      else
         @error "Not implemented for ncells = $ncells."
      end

   # NonAMR Cartesian mesh
   grid1 = range(coordmin[dim_[1]], coordmax[dim_[1]], length=ncells[dim_[1]])
   grid2 = range(coordmin[dim_[2]], coordmax[dim_[2]], length=ncells[dim_[2]])

   ns1 = 14 # number of seeding points on the dayside
   ns2 = 10 # number of seeding points on the nightside
   ns4 = 4  # number of seeding points in the polar region
   seeds = Matrix{Float64}(undef, 2, ns1+2*ns2+2*ns4)

   for i in 1:ns1
      seeds[1,i] = coordmax[dim_[1]] / ns1 * (i - 1)
      seeds[2,i] = 0.0
   end

   for i in 1:ns2
      seeds[1,ns1+i] = -10RE
      seeds[2,ns1+i] = coordmin[dim_[2]] +
         (coordmax[dim_[2]] - coordmin[dim_[2]]) / ns2 * (i - 1)
      seeds[1,ns1+ns2+i] = -30RE
      seeds[2,ns1+ns2+i] = coordmin[dim_[2]] +
         (coordmax[dim_[2]] - coordmin[dim_[2]]) / ns2 * (i - 1)
   end

   for i in 1:ns4
      seeds[1,ns1+2*ns2+i] = -20RE + 10RE*(i - 1)
      seeds[2,ns1+2*ns2+i] = 30RE
      seeds[1,ns1+2*ns2+ns4+i] = -20RE + 10RE*(i - 1)
      seeds[2,ns1+2*ns2+ns4+i] = -30RE
   end

   b = meta["vg_b_vol"]
   b1 = reshape(b[dim_[1],:], ncells[dim_[1]], ncells[dim_[2]])
   b2 = reshape(b[dim_[2],:], ncells[dim_[1]], ncells[dim_[2]])

   for i = axes(seeds,2)
      startx, starty = seeds[:,i]
      x1, y1 = trace(b1, b2, startx, starty, grid1, grid2;
         ds=0.5, maxstep=3000, gridtype="ndgrid")
      x1 ./= RE
      y1 ./= RE
      if length(x1) < 5; continue; end
      line = plot(x1, y1, color="w")
      add_arrow(line[1])
   end

end

main()
```