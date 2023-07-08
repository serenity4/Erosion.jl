using Erosion
using Accessors
using FileIO, ImageIO

!(@isdefined generate_terrain) && includet("../utils.jl")

# terrain = generate_terrain((512, 512))
terrain = generate_terrain((1024, 1024), seed = 4)
# heatmap(terrain)

erosion = HydraulicErosionV1(iterations = 400000, droplet_effect_radius = 0.01, seed = 1, erosion_factor = 0.02, droplet_inertia = 0.6, deposition_factor = 2, droplet_capacity = 2, terrain_size = (2.0, 2.0))

erode_and_save(terrain, erosion; terrain = "examples/terrain_v1.png", eroded = "examples/eroded_v1.png")
