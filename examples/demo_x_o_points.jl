# Finding X-points and O-points in a 2D magnetic reconnection configuration.
# This example assumes X-Z meridional plane.
# Note:
# 1. The input B field domain matters for computing the flux function, but I'm entirely sure
# why there are differences.
# 2. In identifying the X-points and O-points, we currently provide two methods: method 1
# needs to set the relative tolerance, while method 2 does not.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot

function main()
   file = "bulk.0001657.vlsv"
   meta = load(file)

   ndims(meta) != 2 && @error "Flux function only works for 2D simulations!"

   meta.ncells[3] == 1 && @error "equatorial plane, no reconnection!"

   b = meta["vg_b_vol"]
   b = reshape(b, 3, meta.ncells[1], meta.ncells[3])
   dx = [meta.dcoord[1], meta.dcoord[3]]

   x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1]) ./ Vlasiator.RE
   z = LinRange(meta.coordmin[3], meta.coordmax[3], meta.ncells[3]) ./ Vlasiator.RE

   xmin_ = searchsortedfirst(x, 7.0)
   xmax_ = searchsortedlast(x, 9.5)
   zmin_ = searchsortedfirst(z, -5.0)
   zmax_ = searchsortedlast(z, 5.0)

   # meshgrid for plotting
   X = [a for a in x[xmin_:xmax_], _ in z[zmin_:zmax_]]
   Z = [b for _ in x[xmin_:xmax_], b in z[zmin_:zmax_]]

   flux = compute_flux_function(view(b,:,xmin_:xmax_,zmin_:zmax_), dx)

   fig, ax = plt.subplots(subplot_kw=Dict("projection"=>"3d"))

   ax.plot_surface(X, Z, flux;
      cmap=matplotlib.cm.turbo,
      linewidth=0, antialiased=false)

   indices_x, indices_o = find_reconnection_points(flux; retol=1e-4, method=2)

   fig, ax = subplots(figsize=(6,10), constrained_layout=true)

   pcolormesh(meta, "proton/vg_v", ax;
      comp=:z, extent=[5, 10, -7.5, 7.5],
      cmap=matplotlib.cm.RdBu_r)
   s1 = ax.scatter(x[indices_x[1,:].+xmin_.-1], z[indices_x[2,:].+zmin_.-1];
      s=50, marker="x", color="tab:gray")
   s2 = ax.scatter(x[indices_o[1,:].+xmin_.-1], z[indices_o[2,:].+zmin_.-1];
      s=12, marker="o", color="tab:brown")

   ax.legend([s1, s2], ["X-point", "O-point"])
end

main()