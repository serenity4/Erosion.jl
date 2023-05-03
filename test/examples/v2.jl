using Plots

terrain = generate_terrain((512, 512); seed)
water, water_flow, velocity, sediment = initialize_maps(terrain)
model = HydraulicErosionV2(1.0)
erode!(terrain, water, water_flow, velocity, sediment, model)
heatmap(terrain)
heatmap(water)
heatmap(norm.(velocity))
heatmap(sediment)

# See https://github.com/JuliaIO/FFMPEG.jl/issues/53
# @gif for time in time_range(model)
#   erode!(terrain, water, water_flow, velocity, sediment, time, model, rainfall)
#   heatmap(water)
# end
