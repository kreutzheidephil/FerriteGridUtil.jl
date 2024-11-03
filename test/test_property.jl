using FerriteGridUtil
using Ferrite
using Test

@testset "property.jl" begin
    nx = 10
    ny = 10
    nz = 10
    order = 1
    corner1 = [-1/2, -1/2, -1/2]
    corner2 = [1/2, 1/2, 1/2]
    grid = Ferrite.generate_grid(Hexahedron, (nx, ny, nz), Vec{3}((-1/2, -1/2, -1/2)), Vec{3}((1/2, 1/2, 1/2)))
    dim = 3
    ipu = Ferrite.Lagrange{RefHexahedron, order}()^dim
    ipp = Ferrite.Lagrange{RefHexahedron, order}()^(dim-2)
    dh = Ferrite.DofHandler(grid)
    add!(dh, :u, ipu)
    add!(dh, :p, ipp)
    close!(dh)
    dofs_per_node = dim + dim - 2
    dofs_corner1 = FerriteGridUtil.get_dofs_from_coord(dh, corner1, dofs_per_node)
    @test dofs_corner1 == collect(1:dofs_per_node)
end