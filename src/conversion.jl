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
function convert_to_makie_mesh(grid::Grid{dim}; displ::Union{Vector{Vec{dim, T}}, Nothing} = nothing, 
                                                cellset::Union{Vector{<:Ferrite.AbstractCell},OrderedSet{Int},String} = getcells(grid)
                                                disconnectcells::Bool = false) where {dim, T}
    return _convert_to_makie_mesh(grid, displ, disconnectcells, cellset)
end

function _convert_to_makie_mesh(grid::Grid{dim}, displ::Union{Vector{Vec{dim, T}}, disconnectcells::Bool, Nothing}, cellset::OrderedSet{Int}) where {dim, T}
    return _convert_to_makie_mesh(grid, displ, disconnectcells, getcells(grid)[collect(cellset)])
end
function _convert_to_makie_mesh(grid::Grid{dim}, displ::Union{Vector{Vec{dim, T}}, disconnectcells::Bool, Nothing}, cellset::String) where {dim, T}
    return _convert_to_makie_mesh(grid, displ, disconnectcells, getcellset(grid, cellset))
end
function _convert_to_makie_mesh(grid::Grid{dim}, displ::Vector{Vec{dim, T}}, disconnectcells::Bool, cellset::Vector{<:Ferrite.AbstractCell}) where {dim, T}
    @assert length(displ) == getnnodes(grid)
    nodes = [(n.x + u) for (n, u) in zip(grid.nodes, displ)]
    return _convert_to_makie_mesh(nodes, Val(disconnectcells), cellset)
end
function _convert_to_makie_mesh(grid::Grid{dim}, ::Nothing, disconnectcells::Bool, cellset::Vector{<:Ferrite.AbstractCell}) where {dim}
    nodes = [n.x for n in grid.nodes]
    return _convert_to_makie_mesh(nodes, Val(disconnectcells), cellset)
end
function _convert_to_makie_mesh(nodes::Vector{Vec{dim, T}}, ::Val{false}, cells::Vector{C}) where {dim, T, <:Ferrite.AbstractCell}
    nodes = _convert_vec_to_makie.(nodes)
    makiepercell, CM = _get_makie_type_data(C)
    makiecells = Vector{CM}(undef, length(grid) * makieperferrite)
    for (i, cell) in enumerate(cells)
        makiecells[1 + (i-1)*makiepercell : i*makiepercell] .= _convert_cell_to_makie(cell)
    end
    # vcat(_convert_cell_to_makie.(cells )...)
    return Makie.GeometryBasics.Mesh(nodes, cells)
end
function _convert_to_makie_mesh(nodes::Vector{Vec{dim, T}}, ::Val{true}, cells::Vector{<:Ferrite.AbstractCell}) where {dim, T}
    makiepercell, CM = _get_makie_type_data(C)
    nodespercell = nnodes(first(cells))
    nodes = Vector{Makie.GeometryBasics.Point{dim,Float64}}(undef, length(cells)*nodespercell)
    makiecells = Vector{CM}(undef, length(grid) * makieperferrite)
    for (i, cell) in enumerate(cells)
        noderange = 1 + (i-1)*nodespercell : i*nodespercell
        nodes[noderange] .= _convert_vec_to_makie.(nodes[[cell.nodes...]]) # Copy nodes of the current cell into ``nodes``
        disconnectedcell = C(Tuple( n for n in noderange )) # Create a new cell using the copied nodes
        makiecells[1 + (i-1)*makiepercell : i*makiepercell] .= _convert_cell_to_makie(disconnectedcell) # Convert the new cell to Makie
    end
    return Makie.GeometryBasics.Mesh(nodes, cells)
end

_convert_vec_to_makie(v::Vec{dim}) where {dim} = Makie.GeometryBasics.Point{dim,Float64}(v...)

_convert_cell_to_makie(cell::Union{Line,QuadraticLine}) = ( Makie.GeometryBasics.NgonFace{2,Int}(Ferrite.vertices(cell)) ,)
_get_makie_type_data(::Union{Type{Line},Type{QuadraticLine}}) = 1, Makie.GeometryBasics.NgonFace{2,Int}

_convert_cell_to_makie(cell::Union{Triangle,QuadraticTriangle}) = ( Makie.GeometryBasics.NgonFace{3,Int}(Ferrite.vertices(cell)) ,)
_get_makie_type_data(::Union{Type{Triangle},Type{QuadraticTriangle}}) = 1, Makie.GeometryBasics.NgonFace{3,Int}

_convert_cell_to_makie(cell::Union{Quadrilateral,QuadraticQuadrilateral}) = ( Makie.GeometryBasics.NgonFace{4,Int}(Ferrite.vertices(cell)) ,)
_get_makie_type_data(::Union{Type{Quadrilateral},Type{QuadraticQuadrilateral}}) = 1, Makie.GeometryBasics.NgonFace{4,Int}

_convert_cell_to_makie(cell::Union{Tetrahedron,QuadraticTetrahedron}) = Tuple( Makie.GeometryBasics.NgonFace{3,Int}(facet) for facet in Ferrite.facets(cell) )
_get_makie_type_data(::Union{Type{Tetrahedron},Type{QuadraticTetrahedron}}) = 4, Makie.GeometryBasics.NgonFace{3,Int}

_convert_cell_to_makie(cell::Union{Hexahedron,QuadraticHexahedron}) = Tuple( Makie.GeometryBasics.NgonFace{4,Int}(facet) for facet in Ferrite.facets(cell) )
_get_makie_type_data(::Union{Type{Hexahedron},Type{QuadraticHexahedron}}) = 6, Makie.GeometryBasics.NgonFace{4,Int}

function _convert_cells_to_makie(grid::Grid{dim}, cellset::Vector{Int}) where {dim}
    return vcat( _convert_cell_to_makie.(get_cells(grid, cellset))...) # TODO: Optimize: Preallocate final vector and write to it (in a loop?) ?
end
