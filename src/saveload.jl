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
    save(grid::Grid{dim, celltype}, filepath::String) where {dim, celltype}

Saves a grid `grid` in an HDF5 file at `filepath`.
"""
function save(grid::Grid{dim, celltype}, filepath::String) where {dim, celltype}
    f = h5open(filepath, "w")
    _save_nodes!(f, grid)
    _save_facetssets!(f, grid)
    _save_cellsets!(f, grid)
    _save_cells!(f, grid)
    close(f)
    return grid
end

function _save_nodes!(f, grid::Grid{dim, celltype}) where {dim, celltype}
    g = create_group(f, "nodes")
    g["coords"] = collect(n.x.data for n in grid.nodes) 
end

function _save_facetssets!(f, grid::Grid{dim, celltype}) where {dim, celltype}
    g = create_group(f, "face sets")
    for (name, set) in pairs(grid.facetsets)
        g[name] = collect( set )
    end
end

function _save_cellsets!(f, grid::Grid{dim, celltype}) where {dim, celltype}
    g = create_group(f, "cell sets")
    for (name, set) in pairs(grid.cellsets)
        g[name] = collect( set )
    end
end

function _save_cells!(f, grid::Grid{dim, celltype}) where {dim, celltype}
    g = create_group(f, string(celltype))
    g["nodeids"] = collect(cell.nodes for cell in grid.cells)
end

###################################################################################################
###################################################################################################
# TODO: implement, test, document

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


function _load_nodes(f)
    g = f["nodes"]
    coords = read(g, "coords")
    dim = length(coords[1])
    nodes = collect( Node(Tensor{1, dim}(values(c))) for c in coords )
    return nodes
end

function _load_facetsets(f)
    g = f["face sets"]
    if isempty(keys(g))
        facetsets = Dict{String,Set{FaceIndex}}()
    else
        facetsets = Dict(
            collect( begin
                # facetids = collect( FacetIndex(values(f)) for f in read(g, key) )
                facetids = collect( FacetIndex(values(f...)) for f in read(g, key) )
                    key => Set(facetids) end
            for key in keys(g) )
                )
    end
    return facetsets
end

function _load_cellsets(f)
    g = f["cell sets"]
    if isempty(keys(g))
        cellsets = Dict{String,Set{Int}}()
    else
        cellsets = Dict(collect( key => Set(read(g, key)) for key in keys(g) ))
    end
    return cellsets
end

function _load_cells(f)
    cells = nothing
    for cellname in celltypenames
        if haskey(f, cellname)
            g = f[cellname]
            celltype = eval(Symbol(cellname))
            cells = collect( celltype(values(n)) for n in read(g, "nodeids") )
            return cells
        end
    end
    if isnothing(cells)
        throw(error("the file does not contain a group for one of the following cell types $celltypenames"))
    end
    return cells
end

"""
    load(filepath::String)

Returns a grid saved usedin `@refpoint` from an HDF5 file saved at `filepath`.
"""
function load(filepath::String)
    f = h5open(filepath, "r")
    nodes = _load_nodes(f)
    facetsets = _load_facetsets(f)
    cellsets = _load_cellsets(f)
    cells = _load_cells(f)
    close(f)
    return Grid(cells, nodes; facetsets=facetsets, cellsets=cellsets)
end