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
@show typeof(cellset)
@show typeof(allcells)
@show typeof(grid.cells)

# Idea is to get the nodes corresponding to a cell set to plot sub domains. Does the order of 
# nodes matter (when passed to convert_to_makie_mesh())?
function _get_nodes_from_cells(cellset::String)
    return grid.nodes[collect(getcellset(grid, cellset))]
end

function _get_nodes_from_cells(cellset::OrderedSet{Int})
    return grid.nodes[collect(cellset)]
end

function _get_nodes_from_cells(cellset::Vector{<:Ferrite.AbstractCell})
    return missing
end

# convert_to_makie_mesh(grid; cellset = allcells)
# convert_to_makie_mesh(grid; cellset = grid.cells)
# convert_to_makie_mesh(grid)

# displacement = collect(zero(Vec{dim}) for i in 1:getnnodes(grid))
# convert_to_makie_mesh(grid, displacement, grid.cells)