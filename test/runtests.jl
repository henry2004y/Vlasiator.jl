using Vlasiator
using Test

@testset "Vlasiator.jl" begin
   @testset "Reading files" begin
      if Sys.iswindows()
         using ZipFile
         r = ZipFile.Reader("data/bulk_vlsv.zip")
         open(r.files[1].name, "w") do io
            write(io, read(r.files[1], String))
         end
      else
         run(`unzip data/bulk_vlsv.zip`)
      end
      filename = "bulk.0000004.vlsv"
      meta = read_meta(filename)
      # Variable strings reading
      varnames = show_variables(meta)
      @test length(varnames) == 7 && varnames[end] == "vg_rhom"
      # Variable info reading
      varinfo = read_variable_info(meta, "proton/vg_rho")
      @test varinfo.unit == "1/m^3"
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
      @test cellids == collect(4:8)
      # Nearest ID with VDF stored
      @test getNearestCellWithVspace(meta, id) == 8

      # AMR level (this should later be replaced with real AMR data!)
      @test get_max_amr_level(meta) == 0
      @test get_amr_level(meta, id) == 0

      # Compare two VLSV files
      @test compare(filename, filename)

      close(meta.fid) # required for Windows
      rm(filename, force=true)
   end

   @testset "Derived variables" begin

   end

   @testset "Plotting" begin

   end
end
