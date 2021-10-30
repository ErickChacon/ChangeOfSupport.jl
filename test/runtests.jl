using ChangeOfSupport
using Test

@testset "RegularKnots" begin
    knots = RegularKnots(-10, 10, 3, 3, 3)
    knotrange = range(knots)
    @test knotrange == -25:5:25
    knots = RegularKnots(-10, 10, 3, 0, 0)
    knotrange = range(knots)
    @test minimum(knotrange) == -10
    @test maximum(knotrange) == 10
end

# @testset "RegularBsplines" begin
#     b = RegularBsplines(-10, 10, 3, 6)
#     @test range(allknots(b)) == -20:5:25
#     range(internalknots(b)) == -10:5:10
# end

