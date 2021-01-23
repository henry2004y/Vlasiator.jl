using Vlasiator
using Test

@testset "Vlasiator.jl" begin
   run(`tar -zxvf data/bulk_vlsv.tar`)
   filename = "bulk.0000004.vlsv"
   meta = read_meta(filename)
   # ID reading
   cellIDs = read_variable(meta, "CellID")
   IDRef = UInt64[10, 9, 8, 7, 2, 1, 3, 4, 5, 6]
   @test cellIDs == IDRef
   indexRef = [6, 5, 7, 8, 9, 10, 4, 3, 2, 1]
   @test meta.cellIndex == indexRef
   rm(filename)
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
end
