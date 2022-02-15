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
y = LinRange(meta.coordmin[3], meta.coordmax[3], meta.ncells[3]) ./ Vlasiator.RE

xmin_ = 450
xmax_ = 500
ymin_ = 470
ymax_ = 670

# meshgrid for plotting
X = [a for a in x, _ in y]
Y = [b for _ in x, b in y]

flux = compute_flux_function(b, dx, nG)

fig, ax = plt.subplots(subplot_kw=Dict("projection"=>"3d"))

ax.plot_surface(X[xmin_:xmax_,ymin_:ymax_], Y[xmin_:xmax_,ymin_:ymax_],
   flux[xmin_:xmax_,ymin_:ymax_];
   cmap=matplotlib.cm.turbo,
   linewidth=0, antialiased=false)

indices_x, indices_o = find_reconnection_points(flux[xmin_:xmax_,ymin_:ymax_], 1e-1)

fig, ax = subplots(figsize=(6,10), constrained_layout=true)
pcolormesh(meta, "proton/vg_rho", ax; extent=[5, 10, -7.5, 7.5])
s1 = ax.scatter(x[indices_x[1,:].+xmin_.-1], y[indices_x[2,:].+ymin_.-1];
   s=12, marker="x", color="tab:gray")
s2 = ax.scatter(x[indices_o[1,:].+xmin_.-1], y[indices_o[2,:].+ymin_.-1];
   s=12, marker="o", color="tab:brown")

ax.legend([s1, s2], ["X-point", "O-point"])
