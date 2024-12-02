testdir = mktempdir()
@testset "Save and load grid with cells of type $(celltype)" for (celltype, grid) in simplegrids
    filepath = joinpath(testdir, "Test_" * string(celltype) * ".h5")
    @test ! isfile(filepath)
    save(grid, filepath)
    @test isfile(filepath)
    grid₂ = load(filepath)
    @test grid.cells == grid₂.cells
    @test grid.nodes == grid₂.nodes
    @test grid.cellsets  == grid₂.cellsets
    @test grid.facetsets == grid₂.facetsets
    @test grid.nodesets  == grid₂.nodesets
end
rm(testdir; recursive=true)
