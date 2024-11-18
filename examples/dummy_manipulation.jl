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

shiftvalue = Vec{3}(i -> 1.0)#rand(Vec{dim})
scalefactor = rand()
refpoint = Vec{3}((0.5,0.5,0.5))#rand(Vec{dim})

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

function LeviCivita(i)
    return 1
end

function LeviCivita(i, j)
    if (i, j) in [(1, 2)]
        return 1
    elseif (i, j) in [(2, 1)]
        return -1
    else
        return 0
    end
end

function LeviCivita(i, j, k)
    if (i, j, k) in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]
        return 1
    elseif (i, j, k) in [(3, 2, 1), (2, 1, 3), (1, 3, 2)]
        return -1
    else
        return 0
    end
end