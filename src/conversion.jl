"""
    convert_to_makie_mesh()

Return a mesh that can be used for plotting with Makie.jl.
"""
function convert_to_makie_mesh()
    return missing
end
# TODO: implement for different dimensions

#=
function convert_grid_for_Makie(grid::Grid{3})
	nodes = [ Makie.GeometryBasics.Point{3,Float64}(n.x...) for n in grid.nodes ]
	cells = vcat([   [Makie.GeometryBasics.NgonFace{3,Int}(facet) for facet in Ferrite.facets(cell)] for cell in grid.cells   ]...)
    return Makie.GeometryBasics.Mesh(nodes, cells)
end

function convert_grid_for_Makie(grid::Grid{3}, displ::Vector{Vec{3,Float64}})
	nodes = [ Makie.GeometryBasics.Point{3,Float64}( (n.x + u)... ) for (n, u) in zip(grid.nodes, displ) ]
	cells = vcat([   [Makie.GeometryBasics.NgonFace{3,Int}(facet) for facet in Ferrite.facets(cell)] for cell in grid.cells   ]...)
    return Makie.GeometryBasics.Mesh(nodes, cells)
end
=#