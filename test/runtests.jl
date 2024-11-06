using FerriteGridUtil
using Test

@testset "FerriteGridUtil.jl" begin
    @testset "Computing properties" begin
        include("property.jl")
    end
    
    @testset "Manipulating" begin
        include("manipulation.jl")
    end

    @testset "Saving and loading" begin
        include("saveload.jl")
    end

    @testset "Converting the grid" begin
        include("conversion.jl")
    end
end

@info "Tests finished"
