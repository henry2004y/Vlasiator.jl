using Vlasiator
using Test

@testset "Vlasiator.jl" begin
   @testset "Reading files" begin
      if Sys.iswindows()
         using ZipFile
         r = ZipFile.Reader("data/bulk_vlsv.zip")
         for file in r.files
            open(file.name, "w") do io
               write(io, read(file, String))
            end
         end
      else
         run(`unzip data/bulk_vlsv.zip`)
      end

      filenames = ["bulk.0000004.vlsv", "bulk.amr.vlsv"]
      meta = read_meta(filenames[1])
      # Variable strings reading
      varnames = show_variables(meta)
      @test length(varnames) == 7 && varnames[end] == "vg_rhom"
      # Variable info reading
      varinfo = read_variable_info(meta, "proton/vg_rho")
      @test varinfo.unit == "1/m^3"
      # Parameter checking
      @test has_parameter(meta, "dt") == true
      # Parameter reading
      t = read_parameter(meta, "time")
      @test t == 8.0
      # ID reading, unsorted
      cellIDs = read_variable(meta, "CellID", false)
      IDRef = UInt64[10, 9, 8, 7, 2, 1, 3, 4, 5, 6]
      @test cellIDs == IDRef
      indexRef = [6, 5, 7, 8, 9, 10, 4, 3, 2, 1]
      @test meta.cellIndex == indexRef
      # ID finding
      loc = [2.0, 0.0, 0.0]
      id = get_cellid(meta, loc)
      coords = get_cell_coordinates(meta, id)
      @test coords == [2.5, 0.0, 0.0]
      @test read_variable_select(meta, "proton/vg_rho", id)[1][1] ≈ 1.77599 atol=1e-5
      # ID in a line
      point1 = [-2.0, 0.0, 0.0]
      point2 = [2.0, 0.0, 0.0]
      cellids, _, _ = get_cell_in_line(meta, point1, point2)
      @test cellids == collect(4:7)
      # Nearest ID with VDF stored
      @test getNearestCellWithVspace(meta, id) == 8

      # velocity space reading
      vcellids, vcellf = read_velocity_cells(meta, 2; pop="proton")
      V = get_velocity_cell_coordinates(meta, vcellids; pop="proton")
      @test V[:,end] == Float32[2.45, 1.95, 1.95]

      # AMR data reading, dccrg grid
      metaAMR = read_meta(filenames[2])
      maxreflevel = get_max_amr_level(metaAMR)
      sliceoffset = abs(metaAMR.ymin)
      idlist, indexlist = getSliceCellID(metaAMR, sliceoffset, maxreflevel,
         ymin=metaAMR.ymin, ymax=metaAMR.ymax)

      data = read_variable(metaAMR, "proton/vg_rho")
      data = refine_data(metaAMR, idlist, data[indexlist], maxreflevel, :y)
      @test sum(data) ≈ 7.690352275026747e8

      # AMR level
      @test get_max_amr_level(metaAMR) == 2
      @test get_amr_level(metaAMR, idlist[1]) == 1

      # Compare two VLSV files
      @test compare(filenames[1], filenames[1])

      # Explicit IO closure required by Windows
      close(meta.fid)
      close(metaAMR.fid)

      for file in filenames
         rm(file, force=true)
      end
   end

   @testset "Derived variables" begin

   end

   @testset "Plotting" begin

   end
end
