using FerriteGridUtil
using Ferrite

nx = 5
ny = 5
nz = 5
order = 1
# x = Vec(-1/2, -1/2, 0.1)
x = [-1/2, -1/2, 1/2]
dim = 3

corners = Dict([ # Points such that: volume==1, center_of_mass==one(Vec{dim})
    1 => [Vec{1}((0.5,)), Vec{1}((1.5,))],
    2 => [Vec{2}((0.5,0.5)), Vec{2}((1.5,1.5))],
    3 => [Vec{3}((0.5,0.5,0.5)), Vec{3}((1.5,1.5,1.5))]
    ])

grid = Ferrite.generate_grid(Hexahedron, (nx, ny, nz), corners[dim][1], corners[dim][2])

ipu = Lagrange{RefHexahedron, order}()^3
ipp = Lagrange{RefHexahedron, order}()
dh = DofHandler(grid)
add!(dh, :u, ipu)
add!(dh, :p, ipp)
close!(dh)

shiftvalue = rand(Vec{dim})
scalefactor = rand()
refpoint = rand(Vec{dim})

grid_shifted = FerriteGridUtil.shift_by(grid, shiftvalue)
grid_scaled = FerriteGridUtil.scale_relative(grid, scalefactor; refpoint)

for (n, n_shifted) in zip(grid.nodes, grid_shifted.nodes)
    if !(n_shifted.x - n.x  ≈ shiftvalue)
        println("failed shift test")
    end
end

for (n, n_scaled) in zip(grid.nodes, grid_scaled.nodes)
    if !(n_scaled.x ≈ refpoint + scalefactor*(n.x - refpoint))
        println("failed scale test")
    end
end

shiftvaluevec = collect(n_shifted.x - n.x for (n, n_shifted) in zip(grid.nodes, grid_shifted.nodes))