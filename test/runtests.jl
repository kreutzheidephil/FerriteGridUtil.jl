using FerriteGridUtil
using Ferrite
using Test

corners = Dict([ # Points such that: volume==1, center_of_mass==one(Vec{dim})
    1 => [Vec{1}((0.5,)), Vec{1}((1.5,))],
    2 => [Vec{2}((0.5,0.5)), Vec{2}((1.5,1.5))],
    3 => [Vec{3}((0.5,0.5,0.5)), Vec{3}((1.5,1.5,1.5))]
    ])
celltypes = Dict([
    1 => (Line, QuadraticLine),
    2 => (Triangle, QuadraticTriangle, Quadrilateral, QuadraticQuadrilateral),
    3 => (Tetrahedron, Hexahedron) # No grid generator!?!: QuadraticTetrahedron, QuadraticHexahedron
    ])
simplegrids = Dict([
    T => generate_grid(T, Tuple(5 for i in 1:dim), corners[dim]...)
    for dim in 1:3 for T in celltypes[dim]])

@testset "FerriteGridUtil.jl" begin

    mixedgrids = Dict{Int,Vector{Grid}}() # TODO: Also have some tests for mixed grids

    # @testset "Computing properties" begin
    #     include("test_property.jl")
    # end
    
    @testset "Manipulating" begin
        include("test_manipulation.jl")
    end

    # @testset "Saving and loading" begin
    #     include("test_saveload.jl")
    # end

    # @testset "Converting the grid" begin
    #     include("test_conversion.jl")
    # end
end

@info "Tests finished"
