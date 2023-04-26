using Erosion
using Erosion: neighborhood, elliptic_falloff, ErosionMetrics, Cell, GridPosition, corners, bilinear_weights
using Test
using ProceduralNoise
using Plots: heatmap

@testset "Erosion.jl" begin
    @testset "Cell" begin
        position = (20.0, 30.0)
        cell = Cell(position)
        @test cell.bottom_left == GridPosition(position)
        @test cell.bottom_right == GridPosition(position .+ (1, 0))
        @test cell.top_left == GridPosition(position .+ (0, 1))
        @test cell.top_right == GridPosition(position .+ 1)
        weights = bilinear_weights(cell, position)
        @test weights[1] == 1.0
        @test sum(weights) == 1
        @test all(≥(0), weights)
        weights = bilinear_weights(cell, position .+ 0.5)
        @test all(==(0.25), weights)
        weights = bilinear_weights(cell, position .+ (0.13, 0.78))
        @test sum(weights) == 1
        @test all(≥(0), weights)
        @test all(bilinear_weights(cell, corners(cell)[i])[i] == 1.0 for i in 1:4)
    end

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
        metrics = ErosionMetrics(0.23, 0.1, 0.24, 0.465)
        @test isa(repr(metrics), String)
    end
end

resolution = (512, 512)
scale = (2, 2) .^ 4
coords = Tuple(0.0:1/r:(1.0 - 1.0/r) for r in resolution)
grid = collect(Iterators.product(coords...))
perlin = Perlin(scale)
noise = Fractal(perlin, octaves = 8)
terrain = noise.(grid)
terrain = remap.(terrain, Ref((minimum(terrain), maximum(terrain))), Ref((0.0, 1.0)))

heatmap(terrain)

erosion = HydraulicErosion(iterations = 100000, droplet_effect_radius = 0.01, seed = 1, erosion_factor = 0.05, deposition_factor = 10)
# result = erode!(copy(terrain), erosion)
let terrain = copy(terrain); display(@profview(erode!(terrain, erosion))); end
let terrain = copy(terrain); display(erode!(terrain, erosion)); heatmap(terrain); end
display(heatmap(terrain)); let terrain = copy(terrain); display(erode!(terrain, erosion)); heatmap(terrain); end
