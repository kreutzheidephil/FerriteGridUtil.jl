
"""
    get_volume()

Return the volume.
"""
function get_volume()
    return missing
end
# TODO: implementation to compute the volume (area in 2D -> same name?) of the whole grid or just a cellset
# -> Compare computation time of computing the volume by integration or from the node coordinates?
# -> How to treat embedded elements and mixed grids?
# -> Idea: generalize to "get_moment" to also compute 1st and 2nd order moments e.g. integrating "1", "x-x̂", "(x-x̂)⊗(x-x̂)"

#= # Computing the volume by integration

function get_volume(grid::Grid{dim}, cellset::AbstractSet{Int})
    # prepare dh, cv
    V = 0.0
    for cc in CellIterator(dh, cellset)
        Ferrite.reinit!(cv, cc) 
        for qp in 1:getnquadpoints(cv)
            V += getdetJdV(cv, qp)
        end
    end
    return V
end
=#

#= # Computing the volume from the nodal coordinates
function get_volume(nodes, ::Ferrite.Triangle)
    a = norm(nodes[1].x .- nodes[2].x)
    b = norm(nodes[1].x .- nodes[3].x)
    c = norm(nodes[2].x .- nodes[3].x)
    s = (a + b + c) / 2
    return sqrt( s*(s-a)*(s-b)*(s-c) )
end

function get_volume(grid::Grid, i::Integer)
    cell = getcells(grid, i)
    return get_volume(getnodes(grid, [cell.nodes...]), cell)
end

function get_volume(grid::Grid, cellset::AbstractSet{Int})
    vol = 0.0
    for i in cellset
        vol +=  get_volume(grid, i)
    end
    return vol
end
get_volume(grid::Grid, cellset::String) = get_volume(grid, getcellset(grid, cellset))
=#

"""
    get_coordinate_limits()

Return the bounds.
"""
function get_coordinate_limits(grid::Grid{3})
    xyz = [ n.x[i] for n in grid.nodes, i in 1:3 ]
    xmin, xmax, ymin, ymax, zmin, zmax = minimum(xyz[:,1]), maximum(xyz[:,1]), minimum(xyz[:,2]), maximum(xyz[:,2]), minimum(xyz[:,3]), maximum(xyz[:,3])
    return ((xmin, xmax), (ymin, ymax), (zmin, zmax))
end
# TODO: implementation to compute the minimum and maximum coordinates in all directions
# -> Is it more efficient to not collect the coordinates in a matrix, but iterate over grid.nodes and check against current min and max?
# -> Generalize to "dim D"

"""
    get_interface_between_sets()

Return the interface.
"""
function get_interface_between_sets(grid::Grid{dim}, set¹::OrderedSet{Int}, set²::OrderedSet{Int}) where {dim}
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

# add get_dofs_from_nodeid() 

# TODO: implementation to compute set of facets connecting two sets of cells
# -> working for all dimensions like this?
# -> return two facet sets, one from the persepctive of each cell set?

"""
    get_dofs_from_coord(dh::Ferrite.DofHandler, x::Vector, dofs_per_node::Int64; radius=1e-3)

Return the degrees of freedom corresponding to a node at coordinate `x`. 
The node must be within a neighbourhood of radius `radius`. The default `radius` is 1e-4. `dofs_per_node` are the
degrees of freedom per node (i.e. 3 for a 3D displacement problem)
"""
function get_dofs_from_coord(dh::Ferrite.DofHandler, x::Vector, dofs_per_node::Int64; radius=1e-4)
    cells = Ferrite.getcells(dh.grid)
    nodeid = nothing
    # dof_range
    @inline get_coords(n) = Ferrite.get_node_coordinate(dh.grid, n)
    # iterate through the cells nodes and find the first 
    # instance of the node id with coordinate in a ball of radius `radius` around x
    for cellid in eachindex(cells) # !! possiblity for errors if cells is not indexed linearly !!
        nodeid = findfirst(_x -> norm(_x - x) < radius, get_coords.(cells[cellid].nodes))
        if !isnothing(nodeid)
            node_dofs = Ferrite.celldofs(dh, cellid)[dofs_per_node*(nodeid-1)+1:dofs_per_node*nodeid]
            return node_dofs
        end
    end
    throw("No node was found with coordinate $x, try increasing radius")
end