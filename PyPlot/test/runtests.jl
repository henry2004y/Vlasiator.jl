
using Vlasiator, VlasiatorPyPlot, PyPlot, LazyArtifacts
using Suppressor: @capture_err, @suppress_out
using REPL.TerminalMenus
using Test

VlasiatorPyPlot.PyPlot.pygui(false)

function simulate_input(keys...)
   keydict = Dict(:up => "\e[A",
                  :down => "\e[B",
                  :enter => "\r")

   new_stdin = Base.BufferStream()
   for key in keys
      if isa(key, Symbol)
         write(new_stdin, keydict[key])
      else
         write(new_stdin, "$key")
      end
   end
   TerminalMenus.terminal.in_stream = new_stdin

   return
end

@testset "VlasiatorPyPlot" begin
   rootpath = artifact"testdata"

   files = joinpath.(rootpath, ("bulk.1d.vlsv", "bulk.2d.vlsv", "bulk.amr.vlsv"))
   meta1 = load(files[1])
   meta2 = load(files[2])
   meta3 = load(files[3])

   @testset "PyPlot" begin
      # 1D
      meta = meta1
      line = plot(meta, "proton/vg_rho")[1]
      @test line.get_ydata() == meta["proton/vg_rho"]
      line = plot(meta, "proton/vg_v", comp=3)[1]
      @test line.get_ydata() == meta["proton/vg_v"][3,:]
      centers = VlasiatorPyPlot.plotmesh(meta, projection="y")
      points = centers.get_offsets()
      @test size(points) == (10, 2)
      fig = plt.figure()
      ax = fig.add_subplot(projection="3d")
      centers = VlasiatorPyPlot.plotmesh(meta, ax; projection="3d")
      points = centers.get_offsets() # only 2D from 3D coords, might be improved
      @test size(points) == (10, 2)
      plt.clf()

      @test_throws ArgumentError pcolormesh(meta, "proton/vg_rho")

      loc = [2.0, 0.0, 0.0]
      v = VlasiatorPyPlot.vdfslice(meta, loc).get_array()
      @test v[786] == 238.24398578141802
      @test_throws ArgumentError VlasiatorPyPlot.vdfslice(meta, loc, species="helium")
      output = @capture_err begin
         v = VlasiatorPyPlot.vdfslice(meta, loc; slicetype=:bperp).get_array()
      end
      @test v[786] == 4.02741885708042e-10

      output = @capture_err begin
         VlasiatorPyPlot.vdfslice(meta, loc; verbose=true)
      end
      @test startswith(output, "[ Info:")

      # 2D
      meta = meta2
      v = pcolormesh(meta, "proton/vg_rho").get_array()
      @test v[end-2] == 999535.8f0 && length(v) == 6300
      v = pcolormesh(meta, "proton/vg_rho", extent=[0.,1.,0.,2.]).get_array()
      @test v[end] == 0.0 && length(v) == 6
      v = pcolormesh(meta, "proton/vg_rho";
         extent=[-2e7, 2e7, -2e7, 2e7], axisunit=SI, colorscale=Log).get_array()
      @test v[end] == 1.00022675f6 && length(v) == 100
      v = pcolormesh(meta, "fg_b").get_array()
      @test v[1] == 3.0058909f-9
      v = pcolormesh(meta, "proton/vg_v", comp=:x, colorscale=SymLog).get_array()
      @test v[2] == -699935.2f0
      v = contour(meta, "proton/vg_rho").get_array()
      @test length(v) == 8 && v[end] == 5.6e6
      v = contourf(meta, "proton/vg_rho").get_array()
      @test length(v) == 8 && v[end] == 5.6e6
      p = streamplot(meta, "proton/vg_v", comp="xy")
      @test typeof(p) == PyPlot.PyObject
      v = quiver(meta, "proton/vg_v", axisunit=SI, stride=1).get_offsets()
      @test size(v) == (6300, 2)
      # should be tested in 3D, but in principle 2D also works
      pArgs = Vlasiator.set_args(meta, "J", EARTH; normal=:y, origin=0.0)
      v = Vlasiator.prep2dslice(meta, "J", :y, 2, pArgs)
      @test v[62,1] == 9.80949f-12

      # 3D AMR
      meta = meta3
      v = pcolormesh(meta, "proton/vg_rho").get_array()
      @test v[255] == 1.0483886f6 && length(v) == 512
      v = pcolormesh(meta, "proton/vg_v").get_array()
      @test v[255] == 99992.586f0
      v = contour(meta, "proton/vg_rho").get_array()
      @test length(v) == 9 && v[end] == 2.1e6
      pArgs = Vlasiator.set_args(meta, "fg_e", EARTH; normal=:y, origin=0.0)
      v = Vlasiator.prep2dslice(meta, "fg_e", :y, 2, pArgs)
      @test v[16,8] == 0.00032270234f0
   end

   @testset "Plot UI" begin
      meta = meta1
      simulate_input(:enter)
      @suppress_out VlasiatorPyPlot.pui(meta; suppress_output=true)
      @test true # if it reaches here, then pass
      meta = meta2
      simulate_input(fill(:down,5), :enter, :enter)
      @suppress_out VlasiatorPyPlot.pui(meta; suppress_output=true)
      @test true # if it reaches here, then pass
   end
end