@testset "get_moment()" begin
    @test isdefined(FerriteGridUtil, :get_moment)
    # @testset "0ᵗʰ moment for grid with $(T) cells" for (T, grid) in simplegrids
    #     @test isapprox(get_moment(grid, 0), 1.0; rtol=1e-8)
    # end
    # @testset "1ˢᵗ moment for grid with $(T) cells" for (T, grid) in simplegrids
    #     onevec = Vec{Ferrite.getspatialdim(grid)}(i -> 1.0)
    #     @test isapprox(get_moment(grid, 1), onevec; rtol=1e-8)
    # end
    # @testset "2ⁿᵈ moment for grid with $(T) cells" for (T, grid) in simplegrids
    #     MomentTensor = SymmetricTensor{2,Ferrite.getspatialdim(grid)}((i,j) -> i==j ? 13/12 : 1)
    #     @test isapprox(get_moment(grid, 2), MomentTensor; rtol=1e-8)
    # end
end

###################################################################################################
###################################################################################################
# TODO: change 2D grid generation to using two vectors, then write test for limits 
@testset "get_coordinate_limits()" begin
    @test isdefined(FerriteGridUtil, :get_coordinate_limits)
    @testset "$(dim)D grid limits with $(T) cells" for (T, grid) in simplegrids
        test_limits = Tuple( (minimum(corners[Ferrite.getspatialdim(grid)][1][i]), maximum(corners[Ferrite.getspatialdim(grid)][1][i])) for i in 1:3 )
        @test isapprox(get_coordinate_limits(grid), test_limits; rtol=1e-8)
    end
end

###################################################################################################
###################################################################################################

@testset "get_interface_between_sets()" begin
    @test isdefined(FerriteGridUtil, :get_interface_between_sets) 
end

@testset "get_dofs_from_coordinate()" begin
    @test isdefined(FerriteGridUtil, :get_dofs_from_coordinate) 
end