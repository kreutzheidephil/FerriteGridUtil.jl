import FerriteGridUtil
using Ferrite

nx = 10
ny = 10
nz = 10
order = 1
x = [0,0,0]

grid = generate_grid(Hexahedron, (nx, ny, nz), Vec{3}((-1/2, -1/2, -1/2)), Vec{3}((1/2, 1/2, 1/2)))

ip = Lagrange{RefHexahedron, order}()^3
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

FerriteGridUtil.test()
# dof = FerriteGridUtil.get_dofs_from_coord(grid, dh, x)
dof = FerriteGridUtil.get_dofs_from_coord(grid, dh, x)#; fieldname=:u)

# (grid::Ferrite.Grid, dh::Ferrite.DofHandler, x::Vector; fieldname::Symbol = missing, radius::Number=1e-3)