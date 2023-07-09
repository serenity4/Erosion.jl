using Erosion
using Accessors
using FileIO, ImageIO

!(@isdefined generate_terrain) && includet("../utils.jl")
seed = 4

terrain = generate_terrain((512, 512); seed)
model = SimulationBasedErosion(20.0)
model = SimulationBasedErosion(20.0; dissolution_constant = 1.9, deposition_constant = 1.9)
display(heatmap(terrain))
result = erode(terrain, model; progress = true)
(; water, water_flow, velocity, sediment) = result.data
eroded = result.terrain
display(heatmap(eroded))
display(heatmap(eroded .- terrain))
heatmap(terrain)
heatmap(water)
heatmap(norm.(velocity))
heatmap(sediment)

# @gif for time in time_range(model)
#   erode!(terrain, water, water_flow, velocity, sediment, time, model, rainfall)
#   heatmap(water)
# end

terrain = generate_terrain((512, 512); seed)
model = SimulationBasedErosion(20.0; dissolution_constant = 1.9, deposition_constant = 1.9)

erode_and_save(terrain, model; terrain = "examples/terrain_2.png", eroded = "examples/eroded_2.png")
