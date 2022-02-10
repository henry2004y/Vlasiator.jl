using Vlasiator, LaTeXStrings, SHA, LazyArtifacts
using Test

group = get(ENV, "TEST_GROUP", :all) |> Symbol

function nanmaximum(x::AbstractArray{T}) where T<:AbstractFloat
   result = convert(eltype(x), NaN)
   for i in x
      if !isnan(i)
         if isnan(result) || i > result
            result = i
         end
      end
   end
   result
end

@testset "Vlasiator.jl" begin
   rootpath = artifact"testdata"

   files = joinpath.(rootpath, ("bulk.1d.vlsv", "bulk.2d.vlsv", "bulk.amr.vlsv"))
   meta1 = load(files[1])
   meta2 = load(files[2])
   meta3 = load(files[3])

   if group in (:read, :all)
      @testset "Reading files" begin
         @test_throws ArgumentError load("data")
         meta = meta1
         @test ndims(meta) == 1
         @test startswith(repr(meta), "File: bulk.1d.vlsv")
         @test size(meta) == 529201
         # Variable strings reading
         @test meta.variable[end] == "vg_boundarytype"
         # Variable info reading
         varinfo = readvariablemeta(meta, "proton/vg_rho")
         @test startswith(repr(varinfo), "Variable in LaTeX")
         @test varinfo.unit == "1/m^3"
         # Velocity mesh display
         @test startswith(repr(meta.meshes["proton"]), "vblocks")
         # Parameter checking
         @test hasparameter(meta, "dt")
         # Parameter reading
         t = readparameter(meta, "time")
         @test t == 10.0
         @test_throws ArgumentError meta["nonsense"]
         # Do-Block syntax
         t = load(files[1]) do meta
            readparameter(meta, "time")
         end
         @test t == 10.0
         # unsorted ID
         @test readvariable(meta, "CellID", false) == 10:-1:1
         indexRef = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
         @test meta.cellindex == indexRef
         # sorted var by default
         @test meta["vg_boundarytype"] == [4, 4, 1, 1, 1, 1, 1, 1, 3, 3]
         # ID finding
         loc = [2.0, 0.0, 0.0]
         id = getcell(meta, loc)
         coords = getcellcoordinates(meta, id)
         @test coords == (3.0, 0.0, 0.0)
         @test readvariable(meta, "proton/vg_rho", id)[1] == 1.2288102f0
         # ID in a line
         point1 = [-4.0, 0.0, 0.0]
         point2 = [4.0, 0.0, 0.0]
         cellids, _, _ = getcellinline(meta, point1, point2)
         @test cellids == 4:7
         point1 = [-10.1, 0.0, 0.0]
         @test_throws DomainError getcellinline(meta, point1, point2)
         point2 = [10.1, 0.0, 0.0]
         @test_throws DomainError getcellinline(meta, point1, point2)

         # Nearest ID with VDF stored
         @test getnearestcellwithvdf(meta, id) == 5

         # velocity space reading
         vcellids, vcellf = readvcells(meta, 5; species="proton")
         V = getvcellcoordinates(meta, vcellids; species="proton")
         @test V[end] == (2.45f0, 1.95f0, 1.95f0)
         @test_throws ArgumentError readvcells(meta, 2)
         f = Vlasiator.flatten(meta.meshes["proton"], vcellids, vcellf)
         @test f[CartesianIndex(26, 20, 20)] == 85.41775f0
         @test getdensity(meta, f) ≈ 1.8255334f0
         @test getdensity(meta, vcellids, vcellf) ≈ 1.8255334f0
         @test getvelocity(meta, f)[1] ≈ 1.0f0 rtol=3e-3
         @test getvelocity(meta, vcellids, vcellf)[1] ≈ 1.0f0 rtol=3e-3
         @test getpressure(meta, f) ≈ zeros(Float32, 6) atol=1e-16
         @test getmaxwellianity(meta, f) ≈ 5.741325243685855 rtol=1e-4

         # AMR data reading, DCCRG grid
         metaAMR = meta3
         sliceoffset = abs(metaAMR.coordmin[2])
         idlist, indexlist = getslicecell(metaAMR, sliceoffset, 2,
            metaAMR.coordmin[2], metaAMR.coordmax[2])

         data = readvariable(metaAMR, "proton/vg_rho")
         dataslice = refineslice(metaAMR, idlist, data[indexlist], :y)
         @test sum(dataslice) ≈ 7.690352275026747e8
         let err = nothing
            try
               getslicecell(metaAMR, sliceoffset, 1, -2., -1.)
            catch err
            end
            @test err isa Exception
         end
         # AMR ID finding
         loc = [12*Vlasiator.RE, 0.0, 0.0]
         @test getcell(metaAMR, loc) == 0x00000000000002d0

         # AMR level
         @test metaAMR.maxamr == 2
         @test getlevel(metaAMR, idlist[1]) == 1

         # DCCRG utilities
         @test isparent(metaAMR, 1)
         @test !isparent(metaAMR, 1080)
         @test getchildren(metaAMR, 1) == getsiblings(metaAMR, 129)
         @test getparent(metaAMR, 129) == 1
         @test_throws ArgumentError getparent(metaAMR, 5)
         @test_throws ArgumentError getsiblings(metaAMR, 1)

         # FS grid
         data = readvariable(metaAMR, "fg_e")
         ncells, namr = metaAMR.ncells, metaAMR.maxamr
         @test size(data) == (3, ncells[1]*namr^2, ncells[2]*namr^2, ncells[3]*namr^2)
         @test data[:,16,8,8] == [-1.3758785f-7, 3.2213068f-4, -3.1518404f-4]


         # Compare two VLSV files
         @test issame(files[1], files[1])
      end
   end

   if group in (:derive, :all)
      @testset "Derived variables" begin
         meta = meta1
         @test meta["Vmag"] |> sortperm == [3, 4, 1, 2, 8, 9, 10, 5, 7, 6]
         @test meta["Protated"][2,2,3] == 1.7952461f-29
         @test meta["Panisotropy"][4] == 1.0147696f0
         @test meta["Tanisotropy"][4] == 1.0147695404826917
         @test meta["Agyrotropy"][1] |> isnan

         meta = meta2
         @test meta["Bmag"][4] == 3.0052159f-9

         @test meta["Emag"][1,10,99] ≈ 2.6120074f-6

         @test meta["VS"] |> nanmaximum == 1.3726345956957596e6

         @test meta["VA"] |> nanmaximum == 2.3202628822256166e8

         @test meta["MA"][end] == 10.700530839822328

         @test meta["MS"][end] == 16.9375888861409

         @test meta["Vpar"][1] == 698735.3f0

         @test meta["Vperp"][1] == 40982.48f0

         @test meta["T"][1] == 347619.9817319378

         @test meta["Pram"][1] == 8.204415428337215e-10

         @test meta["Pb"][1] == 3.59502523919876e-12

         @test meta["Beta"][1] == 1.3359065984817116

         @test meta["BetaStar"][1] == 229.55170154977864

         @test meta["Poynting"][:,10,10] == [-3.677613f-11, 8.859047f-9, 2.4681486f-9]

         @test meta["IonInertial"][1] == 8.578087716535188e7

         @test readvariable(meta, "Larmor", UInt64[1])[1] == 322324.70603759587

         @test meta["Gyroperiod"][1] == 21.834297799454554

         @test Vlasiator.getdata2d(meta, "Gyroperiod") |> ndims == 2

         @test meta["Gyrofrequency"][1] == 0.04579950356933307

         @test meta["J"][1,1000] == 1.4322914301940016e-15

         @test meta["Omegap"][1] == 209.5467447842415

         @test meta["Plasmaperiod"][1] == 0.0047722048893178645
      end
      @testset "VLSV writing" begin
         meta = meta1
         # Obtain unsorted derived variables, workaround #59
         vmag = readvariable(meta, "Vmag", meta.cellid)
         pa = readvariable(meta, "Panisotropy", meta.cellid)
         vars = Vector{Tuple{VecOrMat, String, VarInfo}}(undef, 0)
         push!(vars, (vmag, "vmag", VarInfo("m/s", L"$\mathrm{m}/mathrm{s}$", L"$V$", "")))
         push!(vars, (pa, "panisotropy", VarInfo("", "", "", "")))

         write_vlsv(files[1], "bulk_new.vlsv", vars)
         sha_str = bytes2hex(open(sha1, "bulk_new.vlsv"))
         @test sha_str == "2b209fb022a2db3b013e01606a60c1fac360a79d"

         rm("bulk_new.vlsv", force=true)
      end
   end

   if group in (:utility, :all)
      @testset "Rotation" begin
         v = fill(1/√3, 3)
         θ = π / 4
         R = Vlasiator.getRotationMatrix(v, θ)
         @test inv(R) ≈ R'

         Rᵀ = Vlasiator.rotateTensorToVectorZ(R, v)
         @test Rᵀ ≈ [1/√2 -1/√2 0.0; 1/√2 1/√2 0.0; 0.0 0.0 1.0]
      end
      @testset "Curvature, Divergence" begin
         let dx = ones(Float32, 3)
            # 2D
            A = ones(Float32, 3,3,1,3)
            @test sum(Vlasiator.curl(dx, A)) == 0.0
            # 3D
            A = ones(Float32, 3,3,3,3)
            @test sum(Vlasiator.curl(dx, A)) == 0.0
            @test sum(Vlasiator.divergence(dx, A)) == 0.0
         end
      end
      @testset "Gradient" begin
         let dx = ones(Float32, 3)
            A = ones(Float32, 3,3,3)
            @test sum(Vlasiator.gradient(dx, A)) == 0.0
         end
      end
   end

   if group in (:vtk, :all)
      @testset "VTK" begin
         meta = meta2 # no amr
         data, ghostType = Vlasiator.fillmesh(meta, ["proton/vg_rho"])
         @test size(data[1][1]) == (1, 63, 100, 1)

         meta = meta3 # amr
         write_vtk(meta, vars=["proton/vg_rho", "fg_b", "proton/vg_v"])
         sha_str = bytes2hex(open(sha1, "bulk.amr_3.vti"))
         @test sha_str == "8a2bb0a15c5dcc329f88821036df840a86eef9d5"

         filesaved = ["bulk.amr.vthb", "bulk.amr_1.vti", "bulk.amr_2.vti", "bulk.amr_3.vti"]
         rm.(filesaved, force=true)
      end
   end

   if group in (:log, :all)
      @testset "Log" begin
         file = joinpath(rootpath, "logfile.txt")
         timestamps, speed = readlog(file)
         @test length(speed) == 50 && speed[end] == 631.2511f0
      end
   end

   if group in (:plot, :all)
      @testset "PyPlot" begin
         using PyPlot
         ENV["MPLBACKEND"]="agg" # no GUI
         # 1D
         meta = meta1
         line = plot(meta, "proton/vg_rho")[1]
         @test line.get_ydata() == meta["proton/vg_rho"]
         line = plot(meta, "proton/vg_v", comp=:3)[1]
         @test line.get_ydata() == meta["proton/vg_v"][3,:]
         centers = plotmesh(meta, projection="y")
         points = centers.get_offsets()
         @test size(points) == (10, 2)
         fig = plt.figure()
         ax = fig.add_subplot(projection="3d")
         centers = plotmesh(meta, projection="3d")
         points = centers.get_offsets() # only 2D from 3D coords, might be improved
         @test size(points) == (10, 2)
         close(fig)

         @test_throws ArgumentError pcolormesh(meta, "proton/vg_rho")

         loc = [2.0, 0.0, 0.0]
         var = vdfslice(meta, loc).get_array()
         @test var[786] == 238.24398578141802
         @test_throws ArgumentError vdfslice(meta, loc, species="helium")

         # 2D
         meta = meta2
         var = pcolormesh(meta, "proton/vg_rho").get_array()
         @test var[end-2] == 999535.8f0 && length(var) == 6300
         var = pcolormesh(meta, "proton/vg_rho", extent=[0,1,0,2]).get_array()
         @test var[end] == 1.0000151f6 && length(var) == 2500
         var = pcolormesh(meta, "proton/vg_rho";
            extent=[-2e7, 2e7, -2e7, 2e7], axisunit=SI).get_array()
         @test var[end] == 1.00022675f6 && length(var) == 100
         var = pcolormesh(meta, "fg_b").get_array()
         @test var[1] == 3.0058909f-9
         var = pcolormesh(meta, "proton/vg_v", comp=:x, colorscale=SymLog).get_array()
         @test var[2] == -699935.2f0
         p = streamplot(meta, "proton/vg_v", comp="xy")
         @test typeof(p) == PyPlot.PyObject
         var = quiver(meta, "proton/vg_v", axisunit=SI, stride=1).get_offsets()
         @test size(var) == (6300, 2)

         # 3D AMR
         meta = meta3
         var = pcolormesh(meta, "proton/vg_rho").get_array()
         @test var[255] == 1.0483886f6 && length(var) == 512
         var = pcolormesh(meta, "proton/vg_v").get_array()
         @test var[255] == 99992.586f0
         @test_throws ArgumentError pcolormesh(meta, "fg_b")
      end

      @testset "Plots" begin
         isdefined(Main, :PyPlot) && import Vlasiator.vdfslice
         include("../src/plot/plots.jl")
         RecipesBase.is_key_supported(k::Symbol) = true
         # 1D
         meta = meta1
         rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), meta, "proton/vg_rho")
         @test getfield(rec[1], 1)[:seriestype] == :line &&
            rec[1].args[1] isa LinRange
#=
         rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(),
            Vlasiator.VDFSlice((meta,[0.0,0.0,0.0])))
         @test getfield(rec[1], 1)[:seriestype] == :histogram2d

         # 2D
         meta = meta2
         rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), meta, "proton/vg_rho")
         @test getfield(rec[1], 1)[:seriestype] == :heatmap &&
            rec[1].args[1] isa LinRange
            =#
      end
   end
   for meta in (meta1, meta2, meta3)
      close(meta.fid) # required for Windows?
   end
end
