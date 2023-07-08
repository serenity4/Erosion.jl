using Erosion
using Accessors
using FileIO, ImageIO

!(@isdefined generate_terrain) && includet("../utils.jl")
seed = 4

using Erosion: initialize_maps

terrain = generate_terrain((512, 512); seed)
water, water_flow, velocity, sediment = initialize_maps(HydraulicErosionV2, terrain)
model = HydraulicErosionV2(20.0)
model = HydraulicErosionV2(20.0; dissolution_constant = 1.9, deposition_constant = 1.9)
display(_heatmap(terrain))
eroded = erode!(copy(terrain), water, water_flow, velocity, sediment, model)
display(_heatmap(eroded))
display(_heatmap(eroded .- terrain))
heatmap(terrain)
heatmap(water)
heatmap(norm.(velocity))
heatmap(sediment)

# @gif for time in time_range(model)
#   erode!(terrain, water, water_flow, velocity, sediment, time, model, rainfall)
#   heatmap(water)
# end

terrain = generate_terrain((512, 512); seed)
model = HydraulicErosionV2(20.0; dissolution_constant = 1.9, deposition_constant = 1.9)

erode_and_save(terrain, model; terrain = "examples/terrain_v2.png", eroded = "examples/eroded_v2.png")
