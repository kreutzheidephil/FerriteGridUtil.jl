# TODO: implement saving and loading grids
# -> to an HDF5 file -> relatively simple, but some work to generalize for all kinds of elements etc.
# -> one of the "NASTRAN® Bulk Data" file types: .nas, .bdf, .nastran, .dat to allow import to COMSOL: https://www.comsol.com/support/learning-center/article/76161
#       |-> Not sure how these files are structured... Gmsh should be able to export to a .bdf file !!! -> A function to add to Gmsh could even be used to export to other types *_*
# -> .stl file for whole grid and only a cellset?

"""
    save()

Saves the grid to a file.
"""
function save()
    return missing    
end

#=
function save_grid(grid::Grid{3,Tetrahedron}, filepath::String)
    f = h5open(filepath, "w")

    g = create_group(f, "nodes")
    g["coords"] = collect(n.x.data for n in grid.nodes)

    g = create_group(f, "tetrahedra")
    g["nodeids"] = collect(cell.nodes for cell in grid.cells)

    g = create_group(f, "face sets")
    for (name, set) in grid.facetsets
        g[name] = collect( f.idx for f in set )
    end

    g = create_group(f, "cell sets")
    for (name, set) in grid.cellsets
        g[name] = collect( set )
    end

    close(f)
    return grid
end
=#

"""
    load()

loads a grid from a file.
"""
function load()
    return missing    
end

#=
function read_grid(filepath::String)
    f = h5open(filepath, "r")

    g = f["nodes"]
    coords = read(g, "coords")
    dim = length(coords[1])
    nodes = collect( Node(Tensor{1,dim}(values(c))) for c in coords )
    
    if haskey(f, "tetrahedra")
        g = f["tetrahedra"]
        cells = collect( Tetrahedron(values(n)) for n in read(g, "nodeids") )
    end

    g = f["face sets"]
    if isempty(keys(g))
        facetsets = Dict{String,Set{FaceIndex}}()
    else
        facetsets = Dict(collect( begin
                facetids = collect( FacetIndex(values(f)) for f in read(g, key) )
                key => Set(facetids) end 
            for key in keys(g) ))
    end

    g = f["cell sets"]
    if isempty(keys(g))
        cellsets = Dict{String,Set{Int}}()
    else
        cellsets = Dict(collect( key => Set(read(g, key)) for key in keys(g) ))
    end

    close(f)
    return Grid(cells, nodes; facetsets=facetsets, cellsets=cellsets)
end
=#