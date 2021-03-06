using Vlasiator, SHA
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
   if Sys.iswindows()
      using ZipFile
      r = ZipFile.Reader("data/testdata.zip")
      for file in r.files
         open(file.name, "w") do io
            write(io, read(file, String))
         end
      end
   else
      run(`unzip data/testdata.zip`)
   end

   filenames = ("bulk.1d.vlsv", "bulk.2d.vlsv", "bulk.amr.vlsv")

   if group in (:read, :all)
      @testset "Reading files" begin
         @test_throws ArgumentError readmeta("data")
         meta = readmeta(filenames[1])
         @test ndims(meta) == 1
         @test startswith(repr(meta), "filename = bulk.1d.vlsv")
         # Variable strings reading
         @test meta.variable[end] == "vg_rhom"
         # Variable info reading
         varinfo = readvariablemeta(meta, "proton/vg_rho")
         @test startswith(repr(varinfo), "var in LaTeX")
         @test varinfo.unit == "1/m^3"
         # Parameter checking
         @test hasparameter(meta, "dt") == true
         # Parameter reading
         t = readparameter(meta, "time")
         @test t == 8.0
         @test_throws ArgumentError readvariable(meta, "nonsense")
         # unsorted ID
         @test readvariable(meta, "CellID", false) == UInt64[10, 9, 8, 7, 2, 1, 3, 4, 5, 6]
         indexRef = [6, 5, 7, 8, 9, 10, 4, 3, 2, 1]
         @test meta.cellIndex == indexRef
         # sorted var by default
         @test readvariable(meta, "vg_boundarytype") == [4, 4, 1, 1, 1, 1, 1, 1, 3, 3]
         # ID finding
         loc = [2.0, 0.0, 0.0]
         id = getcell(meta, loc)
         coords = getcellcoordinates(meta, id)
         @test coords == [2.5, 0.0, 0.0]
         @test readvariable(meta, "proton/vg_rho", id)[1,1] ≈ 1.77599 atol=1e-5
         # ID in a line
         point1 = [-2.0, 0.0, 0.0]
         point2 = [2.0, 0.0, 0.0]
         cellids, _, _ = getcellinline(meta, point1, point2)
         @test cellids == collect(4:7)
         point1 = [-5.1, 0.0, 0.0]
         @test_throws DomainError getcellinline(meta, point1, point2)
         point2 = [5.1, 0.0, 0.0]
         @test_throws DomainError getcellinline(meta, point1, point2)

         @test hasvdf(meta) == true
         # Nearest ID with VDF stored
         @test getnearestcellwithvdf(meta, id) == 8

         # velocity space reading
         vcellids, vcellf = readvcells(meta, 2; pop="proton")
         V = getvcellcoordinates(meta, vcellids; pop="proton")
         @test V[:,end] == Float32[2.45, 1.95, 1.95]
         @test_throws ArgumentError readvcells(meta, 20)

         # AMR data reading, DCCRG grid
         metaAMR = readmeta(filenames[3])
         sliceoffset = abs(metaAMR.coordmin[2])
         idlist, indexlist = getslicecell(metaAMR, sliceoffset;
            ymin=metaAMR.coordmin[2], ymax=metaAMR.coordmax[2])

         data = readvariable(metaAMR, "proton/vg_rho")
         data = refineslice(metaAMR, idlist, data[indexlist], :y)
         @test sum(data) ≈ 7.690352275026747e8
         @test_throws ArgumentError getslicecell(metaAMR, sliceoffset)

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
         @test issame(filenames[1], filenames[1])

         # Explicit IO closure required by Windows
         close(meta.fid)
         close(metaAMR.fid)
      end
   end

   if group in (:derive, :all)
      @testset "Derived variables" begin
         meta = readmeta(filenames[1])
         @test readvariable(meta, "Vmag") |> sortperm == [7, 6, 5, 4, 3, 1, 2, 8, 9, 10]
         close(meta.fid)

         meta = readmeta(filenames[2])
         @test readvariable(meta, "Bmag")[4] ≈ 3.005215661015543e-9

         @test readvariable(meta, "VS") |> nanmaximum == 1.3726345926284448e6

         @test readvariable(meta, "VA") |> nanmaximum == 2.3202627284651753e8

         @test readvariable(meta, "MA")[end] == 10.700530888660992

         @test readvariable(meta, "Vpar")[1] == 698735.3045881701

         @test readvariable(meta, "Vperp")[1] == 40982.48109657114

         @test readvariable(meta, "T")[1] == 347619.97307130496

         @test readvariable(meta, "Beta")[1] == 1.3359065776791028

         Poynting = readvariable(meta, "Poynting")
         @test Poynting[:,10,10] == [-3.677613f-11, 8.859047f-9, 2.4681486f-9]

         @test readvariable(meta, "IonInertial")[1] == 8.584026161906034e7

         @test readvariable(meta, "Larmor")[1] == 142415.61376655987

         #Anisotropy = readvariable(meta, "Anisotropy")

         #Agyrotropy = readvariable(metam "Agyrotropy")

         close(meta.fid)
      end
   end

   if group in (:rotation, :all)
      @testset "Rotation" begin
         using LinearAlgebra
         T = Diagonal([1.0, 2.0, 3.0])
         B = [0.0, 1.0, 0.0]
         Vlasiator.rotateWithB!(T, B)
         @test T == Diagonal([1.0, 3.0, 2.0])
         @test Vlasiator.rotateWithB(T, B) == Diagonal([1.0, 2.0, 3.0])
      end
   end

   if group in (:vtk, :all)
      @testset "VTK" begin
         meta = readmeta(filenames[2]) # no amr
         data, ghostType = Vlasiator.fillmesh(meta, "proton/vg_rho")
         @test size(data[1][1]) == (1, 63, 100, 1)
         close(meta.fid)
         meta = readmeta(filenames[3]) # amr
         write_vtk(meta)
         sha_str = bytes2hex(open(sha1, "bulk.amr_1.vti"))
         @test sha_str == "b127749f30b23d08c814cf169cfaf7fee954bdce"
         close(meta.fid)
         filesaved = ["bulk.amr.vthb", "bulk.amr_1.vti", "bulk.amr_2.vti", "bulk.amr_3.vti"]
         rm.(filesaved, force=true)
      end
   end

   if group in (:log, :all)
      @testset "Log" begin
         file = "logfile.txt"
         timestamps, speed = readlog(file)
         @test length(speed) == 50 && speed[end] == 631.2511f0
         rm.(file, force=true)
      end
   end

   if group in (:plot, :all)
      @testset "Plotting" begin
         using PyPlot
         ENV["MPLBACKEND"]="agg" # no GUI
         # 1D
         meta = readmeta(filenames[1])
         plot(meta, "proton/vg_rho")
         line = gca().lines[1]
         @test line.get_ydata() == readvariable(meta, "proton/vg_rho")
         centers = plotmesh(meta, projection="x")
         points = centers.get_offsets()
         @test size(points) == (1, 2)
         centers = plotmesh(meta, projection="y")
         points = centers.get_offsets()
         @test size(points) == (10, 2)
         centers = plotmesh(meta, projection="z")
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
         p = plot_vdf(meta, loc)
         @test p.get_array()[786] ≈ 229.8948609959216
         @test_throws ArgumentError plot_vdf(meta, loc, pop="helium")
         close(meta.fid)

         # 2D
         meta = readmeta(filenames[2])
         p = pcolormesh(meta, "proton/vg_rho")
         @test p.get_array()[end-2] ≈ 999535.7814279408 && length(p.get_array()) == 6300
         p = pcolormesh(meta, "fg_b")
         @test p.get_array()[1] ≈ 3.0058909f-9
         @test_throws DomainError pcolormesh(meta, "proton/vg_v", op=:x)
         p = streamplot(meta, "proton/vg_v", comp="xy")
         @test typeof(p) == PyPlot.PyObject
         p = quiver(meta, "proton/vg_v", axisunit=SI)
         @test size(p.get_offsets()) == (6300, 2)
         close(meta.fid)

         # 3D AMR
         meta = readmeta(filenames[3])
         p = pcolormesh(meta, "proton/vg_rho")
         @test p.get_array()[255] ≈ 1.04838862e6 && length(p.get_array()) == 512
         @test_throws ArgumentError pcolormesh(meta, "fg_b")
         close(meta.fid)
      end
   end

   for file in filenames
      rm(file, force=true)
   end
end
