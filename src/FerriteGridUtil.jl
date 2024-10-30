module FerriteGridUtil

# Functions to compute properties of a grid
include("property.jl")
export get_volume, get_bounds, get_interface_between_sets

# Functions to manipulate a grid
include("manipulation.jl")
export scale, scale!, shift, shift!

# Functions to save and load a grid
include("save_load.jl")
export save, load

# Functions to convert a grid to a mesh-type of a different package
include("conversion.jl")
export convert_to_makie_mesh

end
