using StaticArrays
using LinearAlgebra
using Erosion: time_range, rainfall_map

function initialize_maps(terrain)
  water = zeros(size(terrain))
  water_flow = zeros(SVector{4,Float64}, size(terrain))
  velocity = zeros(SVector{2,Float64}, size(terrain))
  sediment = zeros(size(terrain))
  water, water_flow, velocity, sediment
end

seed = 1

@testset "Version 2" begin
  model = HydraulicErosionV2(1.0)
  @test model.timestep < 0.5
  terrain = generate_terrain((512, 512); seed)
  water, water_flow, velocity, sediment = initialize_maps(terrain)
  erode!(terrain, water, water_flow, velocity, sediment, model)
  @test !allequal(water) && !allequal(water_flow) && !allequal(velocity)
  @test all(all(x .≥ 0) for x in water_flow)
  @test all(-0.01 .≤ terrain .≤ 1)
  @test all(0 .≤ water .≤ 1)
  @test all(all(.!(isnan.(x))) for x in velocity)
  @test all(0 .≤ sediment .≤ 1)

  model = HydraulicErosionV2(10.0)
  terrain = generate_terrain((512, 512); seed)
  water, water_flow, velocity, sediment = initialize_maps(terrain)
  erode!(terrain, water, water_flow, velocity, sediment, model)
  @test !allequal(water) && !allequal(water_flow) && !allequal(velocity)
  @test all(all(x .≥ 0) for x in water_flow)
  @test all(-0.01 .≤ terrain .≤ 1)
  @test all(-0.05 .≤ water .≤ 1)
  @test all(all(.!(isnan.(x))) for x in velocity)
  @test all(0 .≤ sediment .≤ 1)
end
