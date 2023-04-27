using Erosion
using ProceduralNoise
using Plots: heatmap
using Accessors
using FileIO, ImageIO

!(@isdefined generate_terrain) && includet("utils.jl")

# terrain = generate_terrain((512, 512))
terrain = generate_terrain((1024, 1024))
# heatmap(terrain)

erosion = HydraulicErosion(iterations = 400000, droplet_effect_radius = 0.01, seed = 1, erosion_factor = 0.05, deposition_factor = 10, terrain_size = (2.0, 2.0))

erode_and_save(terrain, erosion)
