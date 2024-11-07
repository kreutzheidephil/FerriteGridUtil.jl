import FerriteGridUtil
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
dofs = FerriteGridUtil.get_dofs_from_coord(dh, x, :u)