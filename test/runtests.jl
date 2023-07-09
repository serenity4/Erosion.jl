using Erosion
using Erosion: neighborhood, elliptic_falloff, ParticleMetrics
using ProceduralNoise
using Test

include("utils.jl")

@testset "Erosion.jl" begin
    include("particle_based.jl")
    include("v2.jl")
end;
