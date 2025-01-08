using FerriteGridUtil
using Ferrite
using OrderedCollections
import GLMakie: Makie, GeometryBasics, Figure, Axis3, Colorbar, Legend, NoShading, wireframe!, mesh!

using Downloads: download
logo_mesh = "logo.geo"
asset_url = "https://raw.githubusercontent.com/Ferrite-FEM/Ferrite.jl/gh-pages/assets/"
isfile(logo_mesh) || download(string(asset_url, logo_mesh), logo_mesh)

grid = togrid(logo_mesh);


# nx = 5
# ny = 5
# nz = 5
# order = 1
# # x = Vec(-1/2, -1/2, 0.1)
# x = [-1/2, -1/2, 1/2]
# dim = 3
# refshape = RefTetrahedron
# shape = Tetrahedron

# grid = Ferrite.generate_grid(shape, (nx, ny, nz))
# addcellset!(grid, "all-cells", x -> true)
# addnodeset!(grid, "rightfacenodes", x -> x[2] <= -1)
# addcellset!(grid, "bottomhalf", x -> x[3] <= 0)

# ipu = Lagrange{refshape, order}()^3
# ipp = Lagrange{refshape, order}()
# dh = DofHandler(grid)
# add!(dh, :u, ipu)
# add!(dh, :p, ipp)
# close!(dh)
# cellset_string = "all-cells"
# cellset_ordereddict = getcellset(grid, "all-cells")
# cellset_vec = getcells(grid)
# # @show typeof(cellset)
# # @show typeof(allcells)
# # @show typeof(cellset_vec)

# u = collect(zero(Vec{dim}) for i in 1:getnnodes(grid))
# for i in getnodeset(grid, "rightfacenodes")
#     u[i] = Vec{dim}([0,rand(-1:0.1:-0.5),0])
# end

@time mesh = convert_to_makie_mesh(grid; disconnectcells=true)#, cellset="bottomhalf")
# mesh = convert_to_makie_mesh(grid; displ=u)#, cellset="bottomhalf")

# fig = Figure()
# ax  = Axis3(fig[1,1]; aspect=:data, title="Grid", xlabel="x [m]", ylabel="y [m]", zlabel="z[m]")
# mesh!(ax, mesh; color=:gray, shading=NoShading)
# wireframe!(ax, mesh; color=:black, linewidth=0.2)
# fig

fig = Makie.Figure()
ax = Makie.Axis(fig[1,1])
Makie.wireframe!(ax, mesh; color=:black)
fig