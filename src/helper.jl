function _get_cellids_by_type(grid::Grid{dim}, cellset::AbstractSet{Int}) where {dim}
    cellidsbytype = Dict([T => OrderedSet{Int}() for T in _get_cell_types_by_dim(Val{dim})])
    for cellid in cellset
        cell = getcells(grid, cellid)
        if ! typeof(cell) in keys(cellidsbytype)
            continue # Exclude not considered cells
        end
        push!(cellidsbytype[typeof(cell)], cellid)
    end
    return cellidsbytype
end
_get_cell_types_by_dim(::Val{1}) = (Line, QuadraticLine)
_get_cell_types_by_dim(::Val{2}) = (Triangle, QuadraticTriangle, Quadrilateral, QuadraticQuadrilateral)
_get_cell_types_by_dim(::Val{3}) = (Tetrahedron, QuadraticTetrahedron, Hexahedron, QuadraticHexahedron, Wedge, Pyramid)
