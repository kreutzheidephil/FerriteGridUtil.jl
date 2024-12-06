# TODO: implement, test, document
# -> return one mesh, or one mesh per type of cell? (One might be interested in treating different types differently, e.g. embedded cells)
#       |-> I am not sure if Makie can treat different elements in the same mesh...
# -> treat quadratic elements differently?
# -> check how it is done in FerriteViz -> will this be needed at all? 
"""
    convert_to_makie_mesh(grid::Grid{dim}; displ::Union{Vector{Vec{dim, T}}, Nothing} = nothing, cellset::Union{Vector{<:Ferrite.AbstractCell},OrderedSet{Int},String} = getcells(grid)) where {dim, T}

Return a mesh that can be used for plotting with Makie.jl. The keyword arguments are:
 - `displ`: a `Vector{Vec{dim, T}}` containing the nodal displacements to plot a deformed grid. The default is undeformed.
 - a `cellset` to plot a subdomain of the grid, the default is the entire `grid`.
"""
function convert_to_makie_mesh(grid::Grid{dim}; displ::Union{Vector{<:Vec{dim, T}}, Nothing} = nothing, cellset::Union{Vector{<:Ferrite.AbstractCell},OrderedSet{Int},String} = getcells(grid)) where {dim, T}
    return _convert_to_makie_mesh(grid, displ, cellset)
end

_convert_to_makie_mesh(grid::Grid{dim}, displ::Union{Vector{<:Vec{dim, T}}, Nothing}, cellset::OrderedSet{Int}) where {dim, T} = _convert_to_makie_mesh(grid, displ, getcells(grid)[collect(cellset)])
_convert_to_makie_mesh(grid::Grid{dim}, displ::Union{Vector{<:Vec{dim, T}}, Nothing}, cellset::String) where {dim, T} = _convert_to_makie_mesh(grid, displ, getcellset(grid, cellset))

function _convert_to_makie_mesh(grid::Grid{dim}, displ::Vector{<:Vec{dim, T}}, cellset::Vector{<:Ferrite.AbstractCell}) where {dim, T}
    @assert length(displ) == getnnodes(grid)
    nodes = [(n.x + u) for (n, u) in zip(grid.nodes, displ)]
    return _convert_to_makie_mesh(nodes, cellset)
end
function _convert_to_makie_mesh(grid::Grid{dim}, ::Nothing, cellset::Vector{<:Ferrite.AbstractCell}) where {dim}
    nodes = [n.x for n in grid.nodes]
    return _convert_to_makie_mesh(nodes, cellset)
end

function _convert_to_makie_mesh(nodes::Vector{<:Vec{dim, T}}, cells::Vector{<:Ferrite.AbstractCell}) where {dim, T}
    nodes = _convert_vec_to_makie.(nodes)
    cells = vcat(_convert_cell_to_makie.(cells )...)
    return Makie.GeometryBasics.Mesh(nodes, cells)
end

_convert_vec_to_makie(v::Vec{dim}) where {dim} = Makie.GeometryBasics.Point{dim,Float64}(v...)
function _convert_cells_to_makie(grid::Grid{dim}, cellset::Vector{Int}) where {dim}
    return vcat( _convert_cell_to_makie.(get_cells(grid, cellset))...) # TODO: Optimize: Preallocate final vector and write to it (in a loop?) ?
end

_convert_cell_to_makie(cell::Union{Line,QuadraticLine}) = [Makie.GeometryBasics.NgonFace{2,Int}(Ferrite.vertices(cell))]
_convert_cell_to_makie(cell::Union{Triangle,QuadraticTriangle}) = [Makie.GeometryBasics.NgonFace{3,Int}(Ferrite.vertices(cell))]
_convert_cell_to_makie(cell::Union{Quadrilateral,QuadraticQuadrilateral}) = [Makie.GeometryBasics.NgonFace{4,Int}(Ferrite.vertices(cell))]
_convert_cell_to_makie(cell::Union{Tetrahedron,QuadraticTetrahedron}) = [Makie.GeometryBasics.NgonFace{3,Int}(facet) for facet in Ferrite.facets(cell)]
_convert_cell_to_makie(cell::Union{Hexahedron,QuadraticHexahedron}) = [Makie.GeometryBasics.NgonFace{4,Int}(facet) for facet in Ferrite.facets(cell)]