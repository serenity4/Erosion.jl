using Erosion
using Erosion: neighborhood, elliptic_falloff, ParticleMetrics, UPLIFT_PRIMITIVE_ASYMMETRIC_RIDGE, ErosionResult
using GeometryExperiments: Point2, Point3, Patch, BezierCurve
using ProceduralNoise
using Test

const P2 = Point2
const P3 = Point3

include("utils.jl")

seed = 1

@testset "Erosion.jl" begin
    include("particle_based.jl")
    include("simulation_based.jl")
    include("uplift.jl")
    include("tectonic_based.jl")
end;
