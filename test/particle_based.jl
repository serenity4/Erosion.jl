@testset "Particle-based erosion" begin
  @testset "Elliptic falloff" begin
      radius = (21.4, 24.9)
      nx, ny = neighborhood(radius)
      @test (nx, ny) == (-22:1:22, -25:1:25)
      radius = (22, 22)
      nx, ny = neighborhood(radius)
      @test (nx, ny) == (-22:1:22, -22:1:22)
      grid = tuple.(nx', ny)
      n = length(grid)
      ws = [elliptic_falloff((0.0, 0.0), point, radius) for point in grid]
      @test count(iszero, ws) < 0.3n
      @test elliptic_falloff((0.0, 0.0), (0.0, 0.0), radius) == 1.0
      @test elliptic_falloff((0.0, 0.0), radius, radius) == 0.0
      @test elliptic_falloff((0.0, 0.0), (0.0, radius[2]), radius) == 0.0
      @test elliptic_falloff((0.0, 0.0), (0.0, radius[2] - 1), radius) > 0.0
      @test elliptic_falloff((0.0, 0.0), (radius[1], 0.0), radius) == 0.0
      @test elliptic_falloff((0.0, 0.0), (radius[1] - 1, 0.0), radius) > 0.0

      radius = (22, 10)
      grid = tuple.(nx', ny)
      n = length(grid)
      ws = [elliptic_falloff((0.0, 0.0), point, radius) for point in grid]
      @test count(iszero, ws) > 0.6n
      @test elliptic_falloff((0.0, 0.0), (0.0, 0.0), radius) == 1.0
      @test elliptic_falloff((0.0, 0.0), radius, radius) == 0.0
      @test elliptic_falloff((0.0, 0.0), (0.0, radius[2]), radius) == 0.0
      @test elliptic_falloff((0.0, 0.0), (0.0, radius[2] - 1), radius) > 0.0
      @test elliptic_falloff((0.0, 0.0), (radius[1], 0.0), radius) == 0.0
      @test elliptic_falloff((0.0, 0.0), (radius[1] - 1, 0.0), radius) > 0.0
  end

  @testset "Erosion metrics" begin
      metrics = ParticleMetrics(0.23, 0.1, 0.24, 0.465)
      @test isa(repr(metrics), String)
  end

  terrain = generate_terrain((512, 512))
  @test all(0 ≤ h ≤ 1 for h in terrain)
  model = ParticleBasedErosion(iterations = 10000)
  result = erode(terrain, model)
  metrics = result.data
  @test metrics.reached_iteration_limit == 0.0
  @test metrics.evaporated < 0.1
  @test metrics.basin > 0.9
  @test metrics.escaped > 0.001
  @test result.terrain ≠ terrain
  @test all(0 ≤ h ≤ 1 for h in result.terrain)
end;
