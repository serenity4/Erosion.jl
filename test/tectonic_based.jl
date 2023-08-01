@testset "Tectonic-based erosion" begin
  terrain = generate_terrain((256, 256); seed)
  model = TectonicBasedErosion(1)
  contour = Patch{BezierCurve,3}(P2[(0.0, 0.0), (0.3, 0.45), (0.3, 0.6), (0.3, 0.7), (0.7, 0.7), (0.8, 0.1), (0.0, 0.0)])
  skeleton = Patch{BezierCurve, 3}(P2[(0.1, 0.1), (0.2, 0.3), (0.3, 0.35), (0.5, 0.4), (0.6, 0.6)])
  data = (; contour, skeleton, uplift = 0.5, radius = 1.5)
  uplift = UpliftPrimitive{Float64}(UPLIFT_PRIMITIVE_ASYMMETRIC_RIDGE, data)
  result = erode(terrain, model, uplift; progress = true)
  @test result isa ErosionResult
  elevation = result.terrain
  (; uplift, drainage_map) = result.data
  @test all(!isnan, uplift)
  @test all(!isnan, drainage)
  @test all(!isnan, elevation)
end;

heatmap(uplift)
heatmap(drainage_map)
heatmap(elevation)
