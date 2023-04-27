module Erosion

using Random: seed!, AbstractRNG, default_rng

abstract type HydraulicErosion end

include("v1.jl")

export HydraulicErosion, erode!, HydraulicErosionV1

end
