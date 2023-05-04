# ---
# title: 1D animation
# id: demo_1d_animation
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to create animation for 1D plots
# ---

This example demonstrates 1D plot animation. `ffmpeg` is required to be installed for saving into mp4.
```julia
using Vlasiator, Glob, PyPlot, Printf

files = glob("bulk.*.vlsv", ".")
# Choose plotting variable
const var = "proton/vg_rho"
const comp = 0 # component of vector, if used

fig = plt.figure(constrained_layout=true)
# Fix axis limits according to data range
ax = plt.axes(xlim=(-10, 10), ylim=(0, 4))

line, = ax.plot([], [], lw=3, marker="*")

function animate(i::Int, files::Vector{String}, var::String, comp::Int, ax, line)
   meta = load(files[i+1])
   x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])
   d = readvariable(meta, var)
   y = ndims(d) == 1 ? d : d[comp,:]
   line.set_data(x, y)

   t = readparameter(meta, "time")
   str_title = @sprintf "t = %4.1fs, var = %s" t var
   ax.set_title(str_title)

   return (line,)
end

# https://matplotlib.org/stable/api/_as_gen/matplotlib.animation.FuncAnimation.html
anim = matplotlib.animation.FuncAnimation(fig, animate, fargs=(files, var, comp, ax, line),
   frames=length(files), blit=true,
   repeat_delay=1000, interval=200)
# Make sure ffmpeg is available!
anim.save("line.mp4", writer="ffmpeg", fps=30)
```