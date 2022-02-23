# Finding X-points and O-points in a 2D magnetic reconnection configuration.
# This example assumes X-Z meridional plane.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot

file = "bulk.0001657.vlsv"

meta = load(file)

nG = 2 # number of ghost cells

ndims(meta) != 2 && @error "Flux function only works for 2D simulations!"

if meta.ncells[3] == 1
   @error "equatorial plane, no reconnection!"
end

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
X = [a for a in x, _ in z]
Z = [b for _ in x, b in z]

flux = compute_flux_function(b, dx, nG)

fig, ax = plt.subplots(subplot_kw=Dict("projection"=>"3d"))

ax.plot_surface(X[xmin_:xmax_,zmin_:zmax_], Z[xmin_:xmax_,zmin_:zmax_],
   flux[xmin_:xmax_,zmin_:zmax_];
   cmap=matplotlib.cm.turbo,
   linewidth=0, antialiased=false)

indices_x, indices_o = find_reconnection_points(flux[xmin_:xmax_,zmin_:zmax_], 5e-3)

fig, ax = subplots(figsize=(6,10), constrained_layout=true)

pcolormesh(meta, "proton/vg_v", ax; comp=:z, extent=[5, 10, -7.5, 7.5])
s1 = ax.scatter(x[indices_x[1,:].+xmin_.-1], z[indices_x[2,:].+zmin_.-1];
   s=50, marker="x", color="tab:gray")
s2 = ax.scatter(x[indices_o[1,:].+xmin_.-1], z[indices_o[2,:].+zmin_.-1];
   s=12, marker="o", color="tab:brown")

ax.legend([s1, s2], ["X-point", "O-point"])
