@testset "scale_relative(), scale_relative!()" begin
    
    # grid_scaled = FerriteGridUtil.scale_relative(grid, scalefactor; refpoint)
    @test isdefined(FerriteGridUtil, :scale_relative) 
    @test isdefined(FerriteGridUtil, :scale_relative!)
    @testset "scale_relative() for grid with $(T) cells" for (T, grid) in simplegrids
        dim = Ferrite.getspatialdim(grid)
        scalefactor = rand(-11:2:11)*rand()
        refpoint = rand(Vec{3})
        grid_scaled = scale_relative(grid, scalefactor; refpoint=refpoint)

        # for (n, n_scaled) in zip(grid.nodes, grid_shifted.nodes)
        #     @test (n_shifted.x - n.x  ≈ shiftvalue)
        # end
        @test get_moment(grid, 0)*scalefactor^dim ≈ get_moment(grid_scaled, 0)
        @test get_moment(grid, 1) ≈ get_moment(grid_scaled, 1) - (refpoint + scalefactor*(onevec{dim} - refpoint))
    end
end

###################################################################################################
###################################################################################################

@testset "shift_by(), shift_by!()" begin
    @test isdefined(FerriteGridUtil, :shift_by) 
    @test isdefined(FerriteGridUtil, :shift_by!)
    
    @testset "shift_by() for grid with $(T) cells" for (T, grid) in simplegrids
        shiftvalue = rand(-11:2:11)*rand(Vec{Ferrite.getspatialdim(grid)})
        grid_shifted = shift_by(grid, shiftvalue)
        # testthis = collect(n_shifted.x - n.x for (n, n_shifted) in zip(grid.nodes, grid_shifted.nodes))
        # @test testthis ≈ shiftvalue
        # TODO: clean up so only one test is performed
        for (n, n_shifted) in zip(grid.nodes, grid_shifted.nodes)
            @test (n_shifted.x - n.x  ≈ shiftvalue)
        end
        @test get_moment(grid, 0) ≈ get_moment(grid_shifted, 0)
        @test get_moment(grid, 1) ≈ get_moment(grid_shifted, 1) - shiftvalue
    end
end