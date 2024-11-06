@testset "get_moment()" begin
    @test isdefined(FerriteGridUtil, :get_moment) 
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
