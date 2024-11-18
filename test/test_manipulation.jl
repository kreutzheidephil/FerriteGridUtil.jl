# @testset "scale_relative(), scale_relative!()" begin
#     # grid_scaled = FerriteGridUtil.scale_relative(grid, scalefactor; refpoint)
#     @test isdefined(FerriteGridUtil, :scale_relative) 
#     @test isdefined(FerriteGridUtil, :scale_relative!)
#     @testset "scale_relative() for grid with $(T) cells" for (T, grid) in simplegrids
#         dim = Ferrite.getspatialdim(grid)
#         scalefactor = rand(1:10)*rand()
#         refpoint = rand(Vec{dim})
#         grid_scaled = scale_relative(grid, scalefactor; refpoint=refpoint)
#         # for (n, n_scaled) in zip(grid.nodes, grid_shifted.nodes)
#         #     @test (n_shifted.x - n.x  ≈ shiftvalue)
#         # end
#         @test get_moment(grid, 0; refpoint)*scalefactor^dim ≈ get_moment(grid_scaled, 0; refpoint)
#         @test get_moment(grid, 1; refpoint)*(scalefactor^(dim+1)) ≈ get_moment(grid_scaled, 1; refpoint)
#         @test get_moment(grid, 2; refpoint)*(scalefactor^(dim+2)) ≈ get_moment(grid_scaled, 2; refpoint)
#     end
# end

###################################################################################################
###################################################################################################

function LeviCivita(i)
    return 1
end

function LeviCivita(i, j)
    if (i, j) in [(1, 2)]
        return 1
    elseif (i, j) in [(2, 1)]
        return -1
    else
        return 0
    end
end

function LeviCivita(i, j, k)
    if (i, j, k) in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]
        return 1
    elseif (i, j, k) in [(3, 2, 1), (2, 1, 3), (1, 3, 2)]
        return -1
    else
        return 0
    end
end

@testset "shift_by(), shift_by!()" begin
    @test isdefined(FerriteGridUtil, :shift_by) 
    @test isdefined(FerriteGridUtil, :shift_by!)
    
    @testset "shift_by() for grid with $(T) cells" for (T, grid) in simplegrids
        dim = Ferrite.getspatialdim(grid)
        shiftvalue = rand(1:10)*rand(Vec{dim})
        grid_shifted = shift_by(grid, shiftvalue)
        ε = Tensor{3, dim}(LeviCivita)
        a = - ε ⋅ shiftvalue
        # @test get_moment(grid, 0) ≈ get_moment(grid_shifted, 0)
        # @test get_moment(grid, 1) ≈ get_moment(grid_shifted, 1) - shiftvalue
        s = get_moment(grid, 1)
        @test get_moment(grid, 2; refpoint=s) + transpose(a) ⋅ a ≈ get_moment(grid_shifted, 2; refpoint=s)
    end
end