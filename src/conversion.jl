# TODO: implement, test, document
# -> return one mesh, or one mesh per type of cell? (One might be interested in treating different types differently, e.g. embedded cells)
#       |-> I am not sure if Makie can treat different elements in the same mesh...
# -> treat quadratic elements differently?
# -> check how it is done in FerriteViz -> will this be needed at all? 
"""
    convert_to_makie_mesh(grid::Grid{dim}, cellset::Union{AbstractSet{Int},String}) where {dim}

Return a mesh that can be used for plotting with Makie.jl.
"""

function convert_to_makie_mesh(grid::Grid{dim}; 
    cellset::Union{Vector{<:Ferrite.AbstractCell},OrderedSet{Int},String} = OrderedSet{Int}(1:getncells(grid))) where {dim}
    println("called 1")
    nodes = [ n.x for n in grid.nodes ]
    @show typeof(cellset)
    @show typeof(_convert_set_to_vec(grid, cellset))
    return convert_to_makie_mesh(nodes, _convert_set_to_vec(grid, cellset)) # calls 2
end

function convert_to_makie_mesh(nodes::Vector{Vec{dim, T}}, cells::Vector{<:Ferrite.AbstractCell}) where {dim, T}
    println("called 2")
    nodes = _convert_vec_to_makie.(nodes)
    cells = vcat(_convert_cell_to_makie.(cells )...)
    return Makie.GeometryBasics.Mesh(nodes, cells)
end

function convert_to_makie_mesh(grid::Grid{dim}, displ::Vector{Vec{dim}}, cellset::Union{Vector{Int},String}) where {dim}
    println("called 3")
    nodes = [ (n.x + u) for (n, u) in zip(grid.nodes, displ) ]
    return convert_to_makie_mesh(nodes, getcells(grid, cellset))
end







_convert_set_to_vec(grid::Grid{dim}, cellset::Vector{<:Ferrite.AbstractCell}) where {dim} = cellset
_convert_set_to_vec(grid::Grid{dim}, cellset::OrderedSet) where {dim} = grid.cells[collect(cellset)]
_convert_set_to_vec(grid::Grid{dim}, cellset::String) where {dim} = _convert_set_to_vec(grid, Ferrite.getcellset(grid,cellset))

_convert_vec_to_makie(v::Vec{dim}) where {dim} = Makie.GeometryBasics.Point{dim,Float64}(v...)
function _convert_cells_to_makie(grid::Grid{dim}, cellset::Vector{Int}) where {dim}
    return vcat( _convert_cell_to_makie.(get_cells(grid, cellset))...)
end
_convert_cell_to_makie(cell::Union{Line,QuadraticLine}) = [Makie.GeometryBasics.NgonFace{2,Int}(Ferrite.vertices(cell))]
_convert_cell_to_makie(cell::Union{Triangle,QuadraticTriangle}) = [Makie.GeometryBasics.NgonFace{3,Int}(Ferrite.vertices(cell))]
_convert_cell_to_makie(cell::Union{Quadrilateral,QuadraticQuadrilateral}) = [Makie.GeometryBasics.NgonFace{4,Int}(Ferrite.vertices(cell))]
_convert_cell_to_makie(cell::Union{Tetrahedron,QuadraticTetrahedron}) = [Makie.GeometryBasics.NgonFace{3,Int}(facet) for facet in Ferrite.facets(cell)]
_convert_cell_to_makie(cell::Union{Hexahedron,QuadraticHexahedron}) = [Makie.GeometryBasics.NgonFace{4,Int}(facet) for facet in Ferrite.facets(cell)]



# function convert_to_makie_mesh(grid::Grid{dim}, cellset::Union{AbstractSet{Int},String}) where {dim}
#     println("called TESTTESTtest")
#     nodes = [ n.x for n in grid.nodes ]
#     println(typeof(nodes))
#     println(typeof(cellset))
#     return convert_to_makie_mesh(nodes, getcells(grid, cellset))
# end
# function convert_to_makie_mesh(grid::Grid{dim}, displ::Vector{Vec{dim}}, cellset::Union{AbstractSet{Int},String}) where {dim}
#     nodes = [ (n.x + u) for (n, u) in zip(grid.nodes, displ) ]
#     return convert_to_makie_mesh(nodes, getcells(grid, cellset))
# end
# function convert_to_makie_mesh(nodes::Vector{Vec{dim, T}}, cells::Vector{<:Ferrite.AbstractCell}) where {dim, T}
#     nodes = _convert_vec_to_makie.( nodes )
#     cells = vcat( _convert_cell_to_makie.( cells )... )
#     return Makie.GeometryBasics.Mesh(nodes, cells)
# end

# _convert_vec_to_makie(v::Vec{dim}) where {dim} = Makie.GeometryBasics.Point{dim,Float64}(v...)
# function _convert_cells_to_makie(grid::Grid{dim}, cellset::AbstractSet{Int}) where {dim}
#     return vcat( _convert_cell_to_makie.( get_cells(grid, cellset) )... )
# end
# _convert_cell_to_makie(cell::Union{Line,QuadraticLine}) = [Makie.GeometryBasics.NgonFace{2,Int}(Ferrite.vertices(cell))]
# _convert_cell_to_makie(cell::Union{Triangle,QuadraticTriangle}) = [Makie.GeometryBasics.NgonFace{3,Int}(Ferrite.vertices(cell))]
# _convert_cell_to_makie(cell::Union{Quadrilateral,QuadraticQuadrilateral}) = [Makie.GeometryBasics.NgonFace{4,Int}(Ferrite.vertices(cell))]
# _convert_cell_to_makie(cell::Union{Tetrahedron,QuadraticTetrahedron}) = [Makie.GeometryBasics.NgonFace{3,Int}(facet) for facet in Ferrite.facets(cell)]
# _convert_cell_to_makie(cell::Union{Hexahedron,QuadraticHexahedron}) = [Makie.GeometryBasics.NgonFace{4,Int}(facet) for facet in Ferrite.facets(cell)]
