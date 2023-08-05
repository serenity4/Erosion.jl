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
          weight = drainage_weight(elevation, point, neighbor, model, (nx, ny), p)
          @test !isnan(weight)
          @test 0 ≤ weight ≤ 1
          weight
        end
        @test sum(weights) < 1.15
      end
    end

    @testset "Convergence of parallel approximation" begin
      model = TectonicBasedErosion(nothing)
      maps = ErosionMaps(elevation, model)
      iterations = 0
      for i in 1:100
        iterations = i
        Threads.@threads for i in 1:nx
          for j in 1:ny
            point = GridPoint(i, j)
            water = compute_drainage(maps.drainage, maps.elevation, point, (nx, ny), model)
            downstream = steepest_neighbor_down(maps.elevation, point, (nx, ny), model)
            downstream_slope = compute_slope(maps.elevation, point, downstream, model, (nx, ny))
            maps.new_drainage[point] = water^0.8 * downstream_slope^2
          end
        end
        isapprox(maps.drainage, maps.new_drainage; rtol = 1e-12) && break
        copyto!(maps.drainage, maps.new_drainage)
      end
      @test 1 < iterations < 10
      low, high = extrema(maps.drainage)
      @test low ≈ 0.0 atol = 1e-14
      @test high < 0.25
    end
  end

  @testset "Erosion algorithm" begin
    contour = Patch{BezierCurve,3}(P2[(0.0, 0.0), (0.3, 0.45), (0.3, 0.6), (0.3, 0.7), (0.7, 0.7), (0.8, 0.1), (0.0, 0.0)])
    skeleton = Patch{BezierCurve, 3}(P2[(0.1, 0.1), (0.2, 0.3), (0.3, 0.35), (0.5, 0.4), (0.6, 0.6)])
    data = (; contour, skeleton, uplift = 0.5, radius = 1.5)
    uplift = UpliftPrimitive{Float64}(UPLIFT_PRIMITIVE_ASYMMETRIC_RIDGE, data)
    model = TectonicBasedErosion(uplift, 5; speed = 5, execution = CPU(parallel = true))
    elevation = generate_terrain((nx, ny); seed)
    result = erode(elevation, model; progress = false)
    @test result isa ErosionResult
    (elevation, (; uplift, drainage)) = (result.terrain, result.data)
    @test all(!isnan, uplift)
    @test all(!isnan, drainage)
    @test all(!isnan, elevation)
    low, high = extrema(elevation)
    @test low ≈ 0.0 atol = 1e-14
    @test high < 1010.0

    @test norm(uplift, Inf) < 1.0
    @test 256*4 ≤ count(≈(Erosion.precipitation(model, (nx, ny))), drainage) ≤ 0.2length(drainage)
    low, high = extrema(drainage)
    @test low ≈ Erosion.precipitation(model, (nx, ny))
    @test high < 115000.0
  end
end;
