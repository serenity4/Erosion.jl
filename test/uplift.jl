@testset "Authoring of uplift maps" begin
  contour = Patch{BezierCurve,3}(P2[(0.0, 0.0), (0.3, 0.45), (0.3, 0.6), (0.3, 0.7), (0.7, 0.7), (0.8, 0.1), (0.0, 0.0)])
  skeleton = Patch{BezierCurve,3}(P2[(0.1, 0.1), (0.2, 0.3), (0.3, 0.35), (0.5, 0.4), (0.6, 0.6)])
  data = (; contour, skeleton, uplift = 0.5, radius = 0.5)
  uplift = UpliftPrimitive{Float64}(UPLIFT_PRIMITIVE_ASYMMETRIC_RIDGE, data)
  @test Erosion.uplift(uplift, P2(0.3, 0.3)) > 0.35
  @test Erosion.uplift(uplift, P2(0.1, 0.1)) ≈ data.uplift
  @test Erosion.uplift(uplift, P2(0.3, 0.35)) ≈ data.uplift
  @test Erosion.uplift(uplift, P2(0.6, 0.6)) ≈ data.uplift
  data = Erosion.materialize_uplift(uplift, (256, 256))
  @test isa(data, Matrix{Float64})

  contour = Patch{BezierCurve,3}(P2f[(0.0, 0.0), (0.3, 0.45), (0.3, 0.6), (0.3, 0.7), (0.7, 0.7), (0.8, 0.1), (0.0, 0.0)])
  skeleton = Patch{BezierCurve,3}(P2f[(0.1, 0.1), (0.2, 0.3), (0.3, 0.35), (0.5, 0.4), (0.6, 0.6)])
  data = (; contour, skeleton, uplift = 0.5f0, radius = 0.5f0)
  uplift = UpliftPrimitive{Float32}(UPLIFT_PRIMITIVE_ASYMMETRIC_RIDGE, data)
  @test isa(Erosion.uplift(uplift, P2f(0.3, 0.3)), Float32)
  @test_broken begin
    data = Erosion.materialize_uplift(uplift, (256, 256))
    @test isa(data, Matrix{Float32})
  end
end;
