
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
The node must be within a neighbourhood of `x` with radius `tol`.
"""
function get_dofs_from_coord(
                            grid::Ferrite.Grid,
                            dh::Ferrite.DofHandler,
                            coord::Vector,
                            fieldname::Symbol = nothing,
                            tol::Number
                            )
        cells = Ferrite.getcells(grid) # all cells from grid are stored in an iterable collection cells[]
        position = nothing
        nodes = Ferrite.getnodes(grid)
        x_to_index = Dict(grid.nodes[n].x => n for n in nodes)
        # iterate through the cells nodes and find the first instance of the node ID corresponding to position coord
        local dofs
        for cellid in 1:length(cells)
            # for nodeid in 1:
            # find cell with node that has correct coordinate, then use celldofs to output the correct dofs
            nodeid = findfirst(_x -> _x == x, cells[cellid].nodes)
            if !isnothing(nodeid)
                dofs = Ferrite.celldofs(dh, i)
                break
            end
        end
        if isanothing(fieldname)
            return dofs
        else
            return missing
        end   
    end

    