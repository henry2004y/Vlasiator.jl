# 1D plot animation.
#
# ffmpeg is required to be installed for saving into mp4.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, Glob, PyPlot

files = glob("bulk.*.vlsv", ".")
# Choose plotting variable
var = "proton/vg_rho"

fig = plt.figure(constrained_layout=true)
# Fix axis limits according to data range
axis = plt.axes(xlim=(-10, 10), ylim=(0, 4))

line, = axis.plot([], [], lw=3, marker="*")

function init()
   line.set_data([], [])
   return (line,)
end

function animate(i::Int, files::Vector{String}, var::String, line)
   meta = load(files[i+1])
   x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])
   y = readvariable(meta, var)
   line.set_data(x, y)

   return (line,)
end

# https://matplotlib.org/stable/api/_as_gen/matplotlib.animation.FuncAnimation.html
anim = matplotlib.animation.FuncAnimation(fig, animate, fargs=(files, var, line),
   init_func=init, frames=length(files), blit=true,
   repeat_delay=1000, interval=200)
# Make sure ffmpeg is available!
anim.save("rho.mp4", writer="ffmpeg", fps=30)