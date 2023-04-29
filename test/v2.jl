using StaticArrays

@testset "Version 2" begin
  terrain = generate_terrain((512, 512))
  water = zeros(size(terrain))
  water_flow = zeros(SVector{4,Float64}, size(terrain))
  velocity = zeros(SVector{2,Float64}, size(terrain))
  sediment = zeros(size(terrain))
  model = HydraulicErosionV2(1.0)
  @test model.timestep < 0.5
  erode!(terrain, water, water_flow, velocity, sediment, model)
  @test !allequal(water) && !allequal(water_flow) && !allequal(velocity)
end
