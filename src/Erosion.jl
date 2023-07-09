module Erosion

using StaticArrays
using GridHelpers: Cell, GridPoint, interpolate_bilinear, bilinear_weights, estimate_gradient, nearest, neighbor, is_outside_grid
using Random: seed!, AbstractRNG, default_rng
using ProceduralNoise: Fractal, Perlin

include("api.jl")

include("particle_based.jl")
include("simulation_based.jl")

export ErosionModel, erode, ParticleBasedErosion, SimulationBasedErosion

end
