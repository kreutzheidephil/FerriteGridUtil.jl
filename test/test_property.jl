@testset "get_moment()" begin
    @test isdefined(FerriteGridUtil, :get_moment)
    @testset "0ᵗʰ moment for grid with $(T) cells" for (T, grid) in simplegrids
        @test isapprox(get_moment(grid, 0), 1.0; rtol=1e-8)
    end
    @testset "1ˢᵗ moment for grid with $(T) cells" for (T, grid) in simplegrids
        onevec = Vec{Ferrite.getspatialdim(grid)}(i -> 1.0)
        @test isapprox(get_moment(grid, 1), onevec; rtol=1e-8)
    end
    @testset "2ⁿᵈ moment for grid with $(T) cells" for (T, grid) in simplegrids
        onevec = Vec{Ferrite.getspatialdim(grid)}(i -> 1.0)
        @test isapprox(get_moment(grid, 2), onevec; rtol=1e-8)
    end
end

###################################################################################################
###################################################################################################

@testset "get_coordinate_limits()" begin
    @test isdefined(FerriteGridUtil, :get_coordinate_limits) 
end

###################################################################################################
###################################################################################################

@testset "get_interface_between_sets()" begin
    @test isdefined(FerriteGridUtil, :get_interface_between_sets) 
end

@testset "get_dofs_from_coord()" begin
    @test isdefined(FerriteGridUtil, :get_dofs_from_coord) 
end