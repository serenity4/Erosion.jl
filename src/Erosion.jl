module Erosion

using StaticArrays
using GridHelpers: Cell, GridPoint, interpolate_bilinear, bilinear_weights, estimate_gradient, nearest, neighbor, is_outside_grid
using Random: seed!, AbstractRNG, default_rng
using ProceduralNoise: Fractal, Perlin

abstract type HydraulicErosion end

include("execution.jl")

include("v1.jl")
include("v2.jl")

export HydraulicErosion, erode!, HydraulicErosionV1, HydraulicErosionV2

end
