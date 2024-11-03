
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
    get_dofs_from_coord()

Return the dofs corresponding to a coordinate. 
The node must be within a neighbourhood of `x` with radius `radius`. The default `radius` is 1e-3 .
"""
function get_dofs_from_coord(grid::Ferrite.Grid, dh::Ferrite.DofHandler, x::Vector; fieldname=nothing, radius=1e-3)
    local dof_position
    if !isnothing(fieldname)
        dof_position = findfirst(name -> name == fieldname, dh.field_names)
        if isnothing(dof_position)
            throw(ArgumentError("Invalid fieldname $fieldname"))
        end
    end
    cells = Ferrite.getcells(grid) # all cells from grid are stored in an iterable collection cells[]
    nodeid = nothing
    local dofs
    @inline get_coords(n) = Ferrite.get_node_coordinate(grid, n)
    # iterate through the cells nodes and find the first 
    # instance of the node ID with coordinate in a ball of radius `tol` around x
    for cellid in eachindex(cells) # !! possiblity for errors if cells is not indexed linearly !!
        nodeid = findfirst(_x -> norm(_x - x) < radius, get_coords.(cells[cellid].nodes))
        if !isnothing(nodeid)
            dofs = Ferrite.celldofs(dh, nodeid)
            break
        end
    end
    if !isnothing(fieldname)
        return dofs[dof_position]
    else
        return dofs
    end
end

function test()
    println("works")
end