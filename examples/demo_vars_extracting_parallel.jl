# Multiple variables alone a line across multiple frames, multi-process version.
#
# To run on a single node,
# julia -p $ncores demo_vars_extracting_parallel_pyplot.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Distributed
@everywhere using Vlasiator, PyPlot, Glob, Printf

@everywhere function init_figure()
   fig, axs = plt.subplots(4, 1; num=myid(),
      figsize=(10, 15), sharex=true, constrained_layout=true)

   axs[end].set_xlim(x1, x2)

   axs[1].set_ylim(ρmin, ρmax)
   axs[2].set_ylim(vmin, vmax)
   axs[3].set_ylim(pmin, pmax)
   axs[4].set_ylim(bmin, bmax)
   for ax in axs
      ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
      ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
      ax.grid(true)
   end

   axs[end].set_xlabel("x [Re]"; fontsize)
   axs[1].set_ylabel("Density [amu/cc]"; fontsize)
   axs[2].set_ylabel("Velocity [km/s]"; fontsize)
   axs[3].set_ylabel("Pressure [nPa]"; fontsize)
   axs[4].set_ylabel("Magnetic field [nT]"; fontsize)
   return fig, axs
end

@everywhere function process(axs, file, cellids)
   isfile("out/"*file[end-8:end-5]*".png") && return

   println("file = $file")
   meta = load(file)

   p_extract = readvariable(meta, "vg_pressure", cellids) .* 1e9 |> vec # [nPa]
   rho_extract = readvariable(meta, "proton/vg_rho", cellids) |> vec
   v_extract = readvariable(meta, "proton/vg_v", cellids)
   vmag2_extract = sum(v_extract.^2, dims=1) |> vec
   pdyn_extract = rho_extract .* Vlasiator.mᵢ .* vmag2_extract .* 1e9 # [nPa]

   bz = readvariable(meta, "vg_b_vol", cellids)[3,:] .* 1e9 #[nT]

   imagnetopause_ = findfirst(<(0.0), bz)

   axs[1].plot(loc, rho_extract ./ 1e6, label="Proton density", color="#1f77b4")
   vl1 = axs[1].vlines(loc[imagnetopause_], ρmin, ρmax;
      colors="r", linestyle="dashed", alpha=0.5)

   axs[2].plot(loc, v_extract[1,:] ./ 1e3, label="Vx", color="#1f77b4")
   vl2 = axs[2].vlines(loc[imagnetopause_], vmin, vmax;
      colors="r", linestyle="dashed", alpha=0.5)

   axs[3].plot(loc, pdyn_extract, label="Dynamic", color="#1f77b4")
   axs[3].plot(loc, p_extract, label="Thermal", color="#ff7f0e")
   vl3 = axs[3].vlines(loc[imagnetopause_], pmin, pmax;
      colors="r", linestyle="dashed", alpha=0.5)

   axs[4].plot(loc, bz, label="Bz", color="#1f77b4")
   hl4 = axs[4].hlines(0.0, loc[1], loc[end], colors="k", linestyle="dashed", alpha=0.2)
   vl4 = axs[4].vlines(loc[imagnetopause_], bmin, bmax;
       colors="r", linestyle="dashed", alpha=0.5)

   str_title = @sprintf "t= %4.1fs" meta.time
   axs[1].set_title(str_title; fontsize)

   axs[2].legend(;loc="upper right", fontsize)
   axs[3].legend(;loc="upper left",  fontsize)
   axs[4].legend(;loc="upper right", fontsize)

   savefig("out/"*file[end-8:end-5]*".png", bbox_inches="tight")

   for ax in axs
      for line in ax.get_lines()
         line.remove()
      end
   end
   for line in (vl1, vl2, vl3, vl4, hl4)
      line.remove()
   end
end

function make_jobs(files)
   for f in files
       put!(jobs, f)
   end
end

@everywhere function do_work(jobs, results)
   fig, axs = init_figure()
   while true
      file = take!(jobs)
      process(axs, file, cellids)
      put!(results, true)
   end
   close(fig)
end

############################################################################################
files = glob("bulk*.vlsv", ".")

const nfile = length(files)

meta = load(files[1])

Re = Vlasiator.Re # Earth radii
const x1, x2 = 8.0, 29.0
point1 = [x1, 0, 0] .* Re
point2 = [x2, 0, 0] .* Re

const cellids, _, _ = getcellinline(meta, point1, point2)

const loc = range(x1, x2, length=length(cellids))

const ρmin, ρmax = 0.0, 10.0     # [amu/cc]
const vmin, vmax = -640.0, 0.0   # [km/s]
const pmin, pmax = 0.0, 1.82     # [nPa]
const bmin, bmax = -25.0, 60.0   # [nT]

const fontsize = 14

@everywhere begin # broadcast parameters
   cellids = $cellids
   x1 = $x1
   x2 = $x2
   ρmin = $ρmin
   ρmax = $ρmax
   vmin = $vmin
   vmax = $vmax
   pmin = $pmin
   pmax = $pmax
   bmin = $bmin
   bmax = $bmax
   fontsize = $fontsize
   loc = $loc
end

const jobs    = RemoteChannel(()->Channel{String}(nfile))
const results = RemoteChannel(()->Channel{Bool}(nfile))

println("Total number of files: $nfile")
println("Running with $(nworkers()) workers...")

@async make_jobs(files) # Feed the jobs channel with all files to process.

@sync for p in workers()
   @async remote_do(do_work, p, jobs, results)
end

n = nfile
@elapsed while n > 0 # wait for all jobs to complete
   take!(results)
   global n = n - 1
end

println("Finished!")
