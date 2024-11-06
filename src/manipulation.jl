# TODO: test, document
"""
    scale_relative()

Return a scaled grid.
"""
function scale_relative(grid::Grid{dim}, scalefactor::Real; refpoint::Vec{dim}=zero(Vec{dim})) where {dim}
    return scale_relative!(deepcopy(grid), scalefactor; refpoint=refpoint)
end

"""
    scale_relative!()

Return a scaled grid.
"""
function scale_relative!(grid::Grid{dim}, scalefactor::Real; refpoint::Vec{dim}=zero(Vec{dim})) where {dim}
    for (i, n) in grid.nodes
        grid.nodes[i] = refpoint + scalefactor*(n.x - refpoint)
    end
    return grid
end

###################################################################################################
###################################################################################################
#TODO: test, document
"""
    shift_by()

Return a shifted grid.
"""
function shift_by(grid::Grid{dim}, shift::Vec{dim}) where {dim}
    return shift_by!(deepcopy(grid), shift)
end

"""
    shift_by!()

Return a shifted grid.
"""
function shift_by!(grid::Grid{dim}, shift::Vec{dim}) where {dim}
    for (i, n) in grid.nodes
        grid.nodes[i] = n.x + shift
    end
    return grid
end

