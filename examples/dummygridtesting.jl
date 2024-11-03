import FerriteGridUtil
using Ferrite

nx = 10
ny = 10
nz = 10
order = 1
x = [-1/2, -1/2, 10.1]

grid = Ferrite.generate_grid(Hexahedron, (nx, ny, nz), Vec{3}((-1/2, -1/2, -1/2)), Vec{3}((1/2, 1/2, 1/2)))

ipu = Lagrange{RefHexahedron, order}()^3
ipp = Lagrange{RefHexahedron, order}()^1
dh = DofHandler(grid)
add!(dh, :u, ipu)
add!(dh, :u, ipp)
close!(dh)
dofs_per_node = 3
dofs = FerriteGridUtil.get_dofs_from_coord(dh, x, dofs_per_node)