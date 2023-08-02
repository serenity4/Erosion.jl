module Erosion

using GeometryExperiments
using StaticArrays
using GridHelpers: Cell, GridPoint, interpolate_bilinear, bilinear_weights, estimate_gradient, nearest, neighbor, neighbors, FourNeighbors, EightNeighbors, is_outside_grid
using Random: seed!, AbstractRNG, default_rng
using ProceduralNoise: Fractal, Perlin
using MLStyle
using LinearAlgebra: norm

include("api.jl")

include("particle_based.jl")
include("simulation_based.jl")
include("uplift.jl")
include("tectonic_based.jl")

export ErosionModel, erode, erode!, ParticleBasedErosion, SimulationBasedErosion, TectonicBasedErosion, ErosionMaps, UpliftPrimitive, uplift

end
