using Erosion
using Erosion: neighborhood, elliptic_falloff, ParticleMetrics
using ProceduralNoise
using Test

include("utils.jl")

@testset "Erosion.jl" begin
    include("particle_based.jl")
    include("simulation_based.jl")
end;
