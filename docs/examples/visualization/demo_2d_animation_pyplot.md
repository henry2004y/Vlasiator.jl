# ---
# title: 2D contour plot animation
# id: demo_2d_contour_animation
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to generate 2D colored contour animation
# ---

This example shows how to plot 2D colored contours with `pcolormesh`.
`ffmpeg` is required to be installed for saving into mp4.
```julia
using Vlasiator, Glob, VlasiatorPyPlot, Printf
using PyCall
@pyimport matplotlib.animation as anim

files = glob("bulk*.vlsv", ".")
const var = "proton/vg_rho"
const comp = 0 # vector component for plotting (if applicable)

fig = plt.figure(figsize=(6.4,5.1), constrained_layout=true)

ax = plt.axes()

const axisunit = EARTH
const colorscale = Log
const addcolorbar = true
# Choose colormap
const cmap = matplotlib.cm.turbo
# Set data plotting range
const vmin = 7.0e4
const vmax = 2.5e6

meta = load(files[1])

pArgs = Vlasiator.set_args(meta, var, axisunit; normal=:none)

norm, ticks = set_colorbar(colorscale, vmin, vmax)

c = let
   x1, x2 = Vlasiator.get_axis(pArgs)
   fakedata = zeros(Float32, length(x2), length(x1))
   ax.pcolormesh(x1, x2, fakedata; norm, cmap)
end

VlasiatorPyPlot.set_plot(c, ax, pArgs, ticks, addcolorbar)

function animate(i::Int, files::Vector{String}, var::String, comp::Int, c)
   meta = load(files[i+1])
   t = readparameter(meta, "time")
   data = Vlasiator.prep2d(meta, var, comp)'
   c.set_array(data)

   str_title = @sprintf "t= %4.1fs" t
   ax.set_title(str_title, fontweight="bold")

   return (c,)
end

# https://matplotlib.org/stable/api/_as_gen/matplotlib.animation.FuncAnimation.html
out = anim.FuncAnimation(fig, animate, fargs=(files, var, comp, c),
   frames=length(files), blit=true,
   repeat_delay=1000, interval=50)
# Make sure ffmpeg is available!
out.save("contour.mp4", writer="ffmpeg", fps=30)
```