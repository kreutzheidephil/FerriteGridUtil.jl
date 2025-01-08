# TODO: implement, test, document
# -> return one mesh, or one mesh per type of cell? (One might be interested in treating different types differently, e.g. embedded cells)
#       |-> I am not sure if Makie can treat different elements in the same mesh...
# -> treat quadratic elements differently?
# -> check how it is done in FerriteViz -> will this be needed at all?


include("cellconversions.jl")

"""
    convert_to_makie_mesh(grid::Grid{dim}; displ::Union{Vector{Vec{dim, T}}, Nothing} = nothing, cellset::Union{Vector{<:Ferrite.AbstractCell},OrderedSet{Int},String} = getcells(grid)) where {dim, T}

Return a mesh that can be used for plotting with Makie.jl. The keyword arguments are:
 - `displ`: a `Vector{Vec{dim, T}}` containing the nodal displacements to plot a deformed grid. The default is undeformed.
 - a `cellset` to plot a subdomain of the grid, the default is the entire `grid`.
"""
function convert_to_makie_mesh(
    grid::Grid{dim}; 
    displ::Union{Vector{Vec{dim, T}}, Nothing} = nothing, 
    cellset::Union{Vector{<:Ferrite.AbstractCell},OrderedSet{Int},String}=getcells(grid),
    disconnectcells::Bool=false
    ) where {dim, T}
    return _convert_to_makie_mesh(grid, displ, disconnectcells, cellset)
end

function _convert_to_makie_mesh(grid::Grid{dim}, displ::Union{Vector{Vec{dim, T}}, Nothing}, disconnectcells::Bool, cellset::OrderedSet{Int}) where {dim, T}
    return _convert_to_makie_mesh(grid, displ, disconnectcells, getcells(grid)[collect(cellset)])
end
function _convert_to_makie_mesh(grid::Grid{dim}, displ::Union{Vector{Vec{dim, T}}, Nothing}, disconnectcells::Bool, cellset::String) where {dim, T}
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
function _convert_to_makie_mesh(nodes::Vector{Vec{dim, T}}, ::Val{false}, cells::Vector{C}) where {dim, T, C<:Ferrite.AbstractCell}
    makienodes = _convert_vec_to_makie.(nodes)
    makiefacepercell, CM = _get_makie_type_data(C)
    makiecells = Vector{CM}(undef, length(cells) * makiefacepercell)
    for (i, cell) in enumerate(cells)
        makiecells[1+(i-1)*makiefacepercell:i*makiefacepercell] .= _convert_cell_to_makie(cell)
    end
    # return vcat(_convert_cell_to_makie.(cells )...)
    return Makie.GeometryBasics.Mesh(makienodes, makiecells)
end

function _convert_to_makie_mesh(nodes::Vector{Vec{dim, T}}, ::Val{true}, cells::Vector{C}) where {dim, T, C<:Ferrite.AbstractCell}
    makiefacepercell, CM = _get_makie_type_data(C)
    nodespercell = Ferrite.nnodes(first(cells))
    makienodes = Vector{Makie.GeometryBasics.Point{dim,Float64}}(undef, length(cells)*nodespercell)
    makiecells = Vector{CM}(undef, length(cells) * makiefacepercell)
    for (i, cell) in enumerate(cells)
        noderange = 1+(i-1)*nodespercell:i*nodespercell
        makienodes[noderange] .= _convert_vec_to_makie.(nodes[[cell.nodes...]]) # Copy nodes of the current cell into ``makienodes``
        disconnected_cell = C(Tuple(n for n in noderange)) # Create a new cell using the copied nodes
        makiecells[1+(i-1)*makiefacepercell:i*makiefacepercell] .= _convert_cell_to_makie(disconnected_cell) # Convert the new cell to Makie
    end
    return Makie.GeometryBasics.Mesh(makienodes, makiecells)
end

# @inline _get_range(quantity::Int64, i::Int64) = 1+(i-1)*quantity:i*quantity
_convert_vec_to_makie(v::Vec{dim}) where {dim} = Makie.GeometryBasics.Point{dim,Float64}(v...)

function _convert_cells_to_makie(grid::Grid{dim}, cellset::Vector{Int}) where {dim}
    return vcat( _convert_cell_to_makie.(get_cells(grid, cellset))...) # TODO: Optimize: Preallocate final vector and write to it (in a loop?) ?
end