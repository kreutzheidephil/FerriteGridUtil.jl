
"""
    get_volume()

Return the volume.
"""
function get_volume()
    return missing
end

"""
    get_volume()

Return the bounds.
"""
function get_bounds()
    return missing
end

"""
    get_volume()

Return the interface.
"""
function get_interface_between_sets()
    return missing
end

"""
    get_dofs_from_coord(dh::Ferrite.DofHandler, x::Vector, dofs_per_node::Int64; radius=1e-3)

Return the degrees of freedom corresponding to a node at coordinate `x`. 
The node must be within a neighbourhood of radius `radius`. The default `radius` is 1e-4. `dofs_per_node` are the
degrees of freedom per node (i.e. 3 for a 3D displacement problem)
"""
function get_dofs_from_coord(dh::Ferrite.DofHandler, x::Vector, dofs_per_node::Int64; radius=1e-4)
    cells = Ferrite.getcells(dh.grid)
    nodeid = nothing
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