using Erosion
using Erosion: neighborhood, elliptic_falloff
using Test

@testset "Erosion.jl" begin
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
end

using ProceduralNoise
using Plots: heatmap

resolution = (512, 512)
scale = (2, 2) .^ 4
coords = Tuple(0.0:1/r:(1.0 - 1.0/r) for r in resolution)
grid = collect(Iterators.product(coords...))
perlin = Perlin(scale)
noise = Fractal(perlin, octaves = 8)
terrain = noise.(grid)

heatmap(terrain)

# result = erode!(copy(terrain), HydraulicErosion())
let terrain = copy(terrain); erode!(terrain, HydraulicErosion()); heatmap(terrain); end
