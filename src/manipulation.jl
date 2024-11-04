"""
    scale_relative()

Return a scaled grid.
"""
function scale_relative(grid::Grid{dim}, scalefactor::Real; xʳᵉᶠ=zero(Vec{dim})) where {dim}
    return scale_relative!(deepcopy(grid), scalefactor; xʳᵉᶠ=xʳᵉᶠ)
end

"""
    scale_relative!()

Return a scaled grid.
"""
function scale_relative!(grid::Grid{dim}, scalefactor::Real; xʳᵉᶠ=zero(Vec{dim})) where {dim}
    return missing
end
# TODO: implementation to scale a grid relative to a point -> shift nodes accordingly

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
    return missing
end
# TODO: implementation to shift a grid by a given vector -> shift nodes accordingly
