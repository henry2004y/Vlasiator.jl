
xloc = range(7.0, 10.0, step=2.35)
yloc = range(-0.0, 5.0, step=2.35)
zloc = 0.0

CIs = CartesianIndices((1:length(xloc), 1:length(yloc)))

fig, axs = plt.subplots(length(xloc), length(yloc), sharex=true, sharey=true)

for i in CIs
   loc = Vlasiator.Re .* [xloc[i[1]], yloc[i[2]], zloc]
   plot_vdf(meta, loc, axs[i]; verbose=false)
end

for a in axs[end,:] a.set_xlabel("vx [km/s]") end
for a in axs[:  ,1] a.set_ylabel("vy [km/s]") end