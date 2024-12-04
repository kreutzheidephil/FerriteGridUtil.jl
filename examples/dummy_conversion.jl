using FerriteGridUtil
using Ferrite, OrderedCollections

nx = 5
ny = 5
nz = 5
order = 1
# x = Vec(-1/2, -1/2, 0.1)
x = [-1/2, -1/2, 1/2]
dim = 3

corners = Dict([ # Points such that: volume==1, center_of_mass==one(Vec{dim})
    1 => [Vec{1}((0.5,)), Vec{1}((1.5,))],
    2 => [Vec{2}((0.5,0.5)), Vec{2}((1.5,1.5))],
    3 => [Vec{3}((0.5,0.5,0.5)), Vec{3}((1.5,1.5,1.5))]
    ])

grid = Ferrite.generate_grid(Hexahedron, (nx, ny, nz), corners[dim][1], corners[dim][2])
addcellset!(grid, "all-cells", x -> true)

ipu = Lagrange{RefHexahedron, order}()^3
ipp = Lagrange{RefHexahedron, order}()
dh = DofHandler(grid)
add!(dh, :u, ipu)
add!(dh, :p, ipp)
close!(dh)
cellset = "all-cells"
allcells = getcellset(grid, "all-cells")
println("type of allcells: $(typeof(allcells))")
println("type of cellset: $(typeof(cellset))")
println("type of grid.cells: $(typeof(grid.cells))")

convert_to_makie_mesh(grid; cellset = cellset)
convert_to_makie_mesh(grid; cellset = allcells)
convert_to_makie_mesh(grid; cellset = grid.cells)
convert_to_makie_mesh(grid)

# function testcall(grid::Grid{dim}, all_cells::OrderedSet) where dim
#     println("successfull call")
#     grid
# end