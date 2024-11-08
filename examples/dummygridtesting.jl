using FerriteGridUtil
using Ferrite

nx = 5
ny = 5
nz = 5
order = 1
# x = Vec(-1/2, -1/2, 0.1)
x = [-1/2, -1/2, 1/2]

grid = Ferrite.generate_grid(Hexahedron, (nx, ny, nz), Vec{3}((-1/2, -1/2, -1/2)), Vec{3}((1/2, 1/2, 1/2)))

ipu = Lagrange{RefHexahedron, order}()^3
ipp = Lagrange{RefHexahedron, order}()
dh = DofHandler(grid)
add!(dh, :u, ipu)
add!(dh, :p, ipp)
close!(dh)
dofs = FerriteGridUtil.get_dofs_from_coordinate(dh, x, :u)
limits=get_coordinate_limits(grid)
corners = Dict([ # Points such that: volume==1, center_of_mass==one(Vec{dim})
    1 => [Vec{1}((0.5,)), Vec{1}((1.5,))],
    2 => [Vec{2}((0.5,0.5)), Vec{2}((1.5,0.5)), Vec{2}((1.5,1.5)), Vec{2}((0.5,1.5))],
    3 => [Vec{3}((0.5,0.5,0.5)), Vec{3}((1.5,1.5,1.5))]
    ])
# Tuple( (minimum(corners[3][][i]), maximum(corners[3][1][i])) for i in 1:3 )
grid = Ferrite.generate_grid(Quadrilateral, (nx, ny), Vec{2}((-1/2, -1/2)), Vec{2}((1/2, 1/2)))
