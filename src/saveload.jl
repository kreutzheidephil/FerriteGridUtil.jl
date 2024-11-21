# TODO: implement, test, document
# -> to an HDF5 file?
#       |-> relatively simple, but some work to generalize for all kinds of elements etc
#       |-> add to FerriteMeshParser.jl?
# -> one of the "NASTRAN® Bulk Data" file types: .nas, .bdf, .nastran, .dat to allow import to COMSOL: https://www.comsol.com/support/learning-center/article/76161
#       |-> Not sure how these files are structured... 
#       |-> Gmsh should be able to export to a .bdf file !!! -> A function to add to Gmsh could even be used to export to other types *_*
#       |-> make this an extension of FerriteGmsh? Or just add it directly to FerriteGmsh? (so far they only translate from Gmsh to Ferrite)
# -> .stl file for whole grid and only a cellset?

"""
    save()

Saves the grid to a file.
"""
function save()
    return missing 
end


function save(grid::Grid{dim, celltype}, filepath::String) where {dim, celltype}
    f = h5open(filepath, "w")

    g = create_group(f, "nodes")
    g["coords"] = collect(n.x.data for n in grid.nodes)

    g = create_group(f, string(celltype))
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


###################################################################################################
###################################################################################################
# TODO: implement, test, document
"""
    load()

loads a grid from a file.
"""
function load()
    return missing    
end

# celltypenames = Dict([
#     "Line" => Line, 
#     "QuadraticLine" => QuadraticLine,
#     "Triangle" => Triangle, 
#     "QuadraticTriangle" => QuadraticTriangle, 
#     "Quadrilateral" => Quadrilateral, 
#     "QuadraticQuadrilateral" => QuadraticQuadrilateral,
#     "Tetrahedron" => Tetrahedron,
#     "Hexahedron" => Hexahedron
#     ])
celltypenames = (
    "Line", 
    "QuadraticLine",
    "Triangle", 
    "QuadraticTriangle", 
    "Quadrilateral", 
    "QuadraticQuadrilateral" ,
    "Tetrahedron",
    "Hexahedron"
    )

function load(filepath::String)
    f = h5open(filepath, "r")

    g = f["nodes"]
    coords = read(g, "coords")
    dim = length(coords[1])
    nodes = collect( Node(Tensor{1,dim}(values(c))) for c in coords )

###########################################
    # loop through cell types
    for cellname in celltypenames 
        if haskey(f, cellname)
            g = f[cellname]
            celltype = eval(Symbol(cellname))
            cells = collect( celltype(values(n)) for n in read(g, "nodeids") )
            break
        end
    end
################################
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

# function _extract_cells()