# TODO: test, document
"""
    get_moment()

Return the volume.
"""
function get_moment(grid::Grid{dim}, order::Int; 
        cellset::Union{String,AbstractSet{Int}}=OrderedSet{Int}(1:getncells(grid)), 
        refpoint::Vec{dim}=zero(Vec{dim})) where {dim}
    if cellset isa String
        cellset = getcellset(grid, cellset)
    end
    cellidsbytype = _get_cellids_by_type(grid, cellset)
    return _compute_moment(grid, Val(order), cellidsbytype, refpoint)
end

function _compute_moment(grid::Grid{dim}, order::Val{ord}, cellidsbytype::Dict{DataType,OrderedSet{Int}}, refpoint::Vec{dim}) where {dim, ord}
    m = _init_moment(order, dim)
    for (T, cellids) in cellidsbytype
        cv = _init_cv(T, order)
        for cellid in cellids
            coords = getcoordinates(grid, cellid)
            reinit!(cv, coords)
            for qp in 1:getnquadpoints(cv)
                dΩ = getdetJdV(cv, qp)
                x  = spatial_coordinate(cv, qp, coords)
                m += _compute_point_moment(x, refpoint, order) * dΩ
            end
        end
    end
    return m
end

_init_moment(::Val{0}, dim::Int) = 0.0
_init_moment(::Val{1}, dim::Int) = zero(Vec{dim})
_init_moment(::Val{2}, dim::Int) = zero(Tensor{2,dim})

function _init_cv(T::Type{<:Ferrite.AbstractCell{RefShape}}, ::Val{ord}) where {RefShape, ord}
    qr  = QuadratureRule{RefShape}(1+ord)
    ipf = Lagrange{RefShape,1}()
    ipg = geometric_interpolation(T)
    return CellValues(qr, ipf, ipg)
end

_compute_point_moment(x::Vec{dim}, x̂::Vec{dim}, ::Val{0}) where {dim} = 1.0
_compute_point_moment(x::Vec{dim}, x̂::Vec{dim}, ::Val{1}) where {dim} = x - x̂
_compute_point_moment(x::Vec{dim}, x̂::Vec{dim}, ::Val{2}) where {dim} = (x - x̂) ⊗ (x - x̂)

###################################################################################################
###################################################################################################
# TODO: test, document
# -> Is it more efficient to not collect the coordinates in a matrix, but iterate over grid.nodes and check against current min and max?
"""
    get_coordinate_limits()

Return the bounds.
"""
function get_coordinate_limits(grid::Grid{dim}) where {dim}
    coords = [ n.x[i] for n in grid.nodes, i in 1:dim ]
    return Tuple( (minimum(coords[:,i]), maximum(coords[:,i])) for i in 1:dim )
end

###################################################################################################
###################################################################################################
# TODO: test, document
# -> working for all dimensions like this?
# -> return two facet sets, one from the perspective of each cell set?
"""
    get_interface_between_sets()

Return the interface.
"""
function get_interface_between_sets(grid::Grid{dim}, set¹::AbstractSet{Int}, set²::AbstractSet{Int}) where {dim}
    top = ExclusiveTopology(grid)
    Γⁱⁿᵗ = OrderedSet{FacetIndex}()
    for cellid¹ in set¹
        cell¹ = grid.cells[cellid¹]
        for facet in 1:nfacets(cell¹)
            facetindex = FacetIndex(cellid¹, facet)
            neighbor = getneighborhood(top, grid, facetindex)
            length(neighbor) == 0 ? continue : nothing
            if neighbor[1][1] in set²
                push!(Γⁱⁿᵗ, facetindex)
            end
        end
    end
    return Γⁱⁿᵗ
end
function get_interface_between_sets(grid, set¹::String, set²::String)
    return get_interface_between_sets(grid, getcellset(grid, set¹), getcellset(grid, set²))
end
function get_interface_between_sets(grid, set¹::String, set²::AbstractSet{Int})
    return get_interface_between_sets(grid, getcellset(grid, set¹), set²)
end
function get_interface_between_sets(grid, set¹::AbstractSet{Int}, set²::String)
    return get_interface_between_sets(grid, set¹, getcellset(grid, set²))
end