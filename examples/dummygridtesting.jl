
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
The node must be within a neighbourhood with radius `tol`.
"""
function get_dofs_from_coord(
                            grid::G,
                            dh::DH,
                            coord::Vector,
                            fieldname::Symbol = nothing,
                            nodecoord2index::Dict,
                            DOF::Int=0
                            ) where {G <: Ferrite.AbstractGrid, DH <: Ferrite.AbstractDofHandler}
        cells = Ferrite.getcells(grid) # all cells from grid are stored in an iterable collection cells[]
        nodeid = nothing
        i = 0
        # iterate through the cells nodes and find the first instance of the node ID corresponding to position coord
        while (isnothing(nodeid) && i < length(cells)) 
            i+=1
            nodeid = findfirst(x -> x == nodecoord2index[coord], cells[i].nodes)
        end
            dofs = celldofs(dh, i) # get the dofs for the cell with node ID i
        if DOF == 3
            return dofs[3*nodeid]
        elseif DOF == 2
            return dofs[3*nodeid-1]
        elseif DOF == 1
            return dofs[3*nodeid-2]
        elseif DOF == 0
            return dofs[3*nodeid-2:3*nodeid]
        else
            error("coord2dof is only for 3D (displacement) problems, the last argument must be either 0, 1, 2 or 3. Here 0 will return all three DOFs (x, y, z). The default value is 0")
        end
        return missing
    end
    