using Erosion
using Erosion: neighborhood, elliptic_falloff, ParticleMetrics, UPLIFT_PRIMITIVE_ASYMMETRIC_RIDGE, ErosionResult, drainage_weight, uplift_map, compute_drainage, compute_slope, steepest_neighbor_down
using GeometryExperiments: Point2, Point3, Point2f, Point3f, Patch, BezierCurve
using GridHelpers: GridPoint, neighbors, EightNeighbors, is_outside_grid
using LinearAlgebra: norm
using ProceduralNoise
using Test
using FileIO, ImageIO
using ImageTransformations

const P2 = Point2
const P3 = Point3
const P2f = Point2f
const P3f = Point3f

include("utils.jl")

seed = 1

@testset "Erosion.jl" begin
    include("particle_based.jl")
    include("simulation_based.jl")
    include("uplift.jl")
    include("tectonic_based.jl")
end;
