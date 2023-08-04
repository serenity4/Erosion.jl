@testset "Tectonic-based erosion" begin
  nx, ny = (256, 256)
  elevation = generate_terrain((nx, ny); seed)

  @testset "Drainage" begin
    @testset "Drainage weights" begin
      p = 4
      for point in GridPoint.([1, 100, 256], [4, 53, 256])
        weights = map(neighbors(point, EightNeighbors())) do neighbor
          is_outside_grid(neighbor, (nx, ny)) && return 0.0
          elevation[neighbor] > elevation[point] || return 0.0
          weight = drainage_weight(elevation, point, neighbor, (nx, ny), p)
          @test !isnan(weight)
          @test 0 ≤ weight ≤ 1
          weight
        end
        @test sum(weights) < 1.15
      end
    end

    model = TectonicBasedErosion(nothing)
    drainage = zeros((nx, ny))
    new_drainage = zeros((nx, ny))
    for i in 1:100
      Threads.@threads for i in 1:nx
        for j in 1:ny
          point = GridPoint(i, j)
          new_drainage[point] = compute_drainage(drainage, elevation, point, (nx, ny), model)
        end
      end
      copyto!(drainage, new_drainage)
      # i % 10 == 1 && display(heatmap(drainage; colormap = :grays))
    end
  end

  contour = Patch{BezierCurve,3}(P2[(0.0, 0.0), (0.3, 0.45), (0.3, 0.6), (0.3, 0.7), (0.7, 0.7), (0.8, 0.1), (0.0, 0.0)])
  skeleton = Patch{BezierCurve, 3}(P2[(0.1, 0.1), (0.2, 0.3), (0.3, 0.35), (0.5, 0.4), (0.6, 0.6)])
  data = (; contour, skeleton, uplift = 0.5, radius = 1.5)
  uplift = UpliftPrimitive{Float64}(UPLIFT_PRIMITIVE_ASYMMETRIC_RIDGE, data)
  model = TectonicBasedErosion(uplift, 5; speed = 5)
  elevation = generate_terrain((nx, ny))
  result = erode(elevation, model; progress = true)
  @test result isa ErosionResult
  elevation = result.terrain
  (; uplift, drainage) = result.data
  @test all(!isnan, uplift)
  @test all(!isnan, drainage)
  @test all(!isnan, elevation)
  @test norm(uplift, Inf) < 1.0
  @test norm(drainage, Inf) < 10.0
  @test norm(elevation, Inf) < 1.0
  @test all(-0.1 .≤ uplift .≤ 1.1)
  @test all(-0.1 .≤ drainage .≤ 1.1)
  @test all(-0.1 .≤ elevation .≤ 1.1)
end;

# heatmap(uplift)
# heatmap(drainage)
# heatmap(elevation)
