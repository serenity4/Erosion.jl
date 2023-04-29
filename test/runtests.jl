using Erosion
using Erosion: neighborhood, elliptic_falloff, ErosionMetricsV1
using ProceduralNoise
using Test

include("utils.jl")

@testset "Erosion.jl" begin
    include("v1.jl")
    include("v2.jl")
end;
