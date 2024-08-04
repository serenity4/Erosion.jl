using ProceduralNoise
using ProceduralNoise: remap
using Random: default_rng, seed!
using ColorTypes

load_uplift(file) = (x -> convert(Float64, convert(Gray, x))).(load(joinpath(@__DIR__, "uplift", file))')
load_elevation(file) = (x -> convert(Float64, convert(Gray, x))).(load(joinpath(@__DIR__, "elevation", file))')
save_uplift(uplift, file) = save(joinpath(@__DIR__, "uplift", file), remap.(uplift, extrema(uplift)..., 0, 1))
save_elevation(elevation, file) = save(joinpath(@__DIR__, "elevation", file), remap.(elevation, extrema(elevation)..., 0, 1))

generate_terrain(resolution; seed = rand(UInt64)) = generate_terrain(Float64, resolution; seed)
function generate_terrain(::Type{T}, resolution; seed = rand(UInt64)) where {T}
    seed!(default_rng(), seed)
    noise = Fractal{Perlin{T}}((2, 2), octaves = 8)
    noise(resolution)
end
