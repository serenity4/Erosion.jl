# Assymetric ridge uplift primitive.

contour = Patch{BezierCurve,3}(P2[(0.0, 0.0), (0.3, 0.45), (0.3, 0.6), (0.3, 0.7), (0.7, 0.7), (0.8, 0.1), (0.0, 0.0)])
skeleton = Patch{BezierCurve, 3}(P2[(0.1, 0.1), (0.2, 0.3), (0.3, 0.35), (0.5, 0.4), (0.6, 0.6)])
data = (; contour, skeleton, uplift = 0.5, radius = 1.5)
p = UpliftPrimitive{Float64}(UPLIFT_PRIMITIVE_ASYMMETRIC_RIDGE, data)

function uplift!(axis, fig, p::UpliftPrimitive)
  n = 512
  xs = ys = 0:(1/(n - 1)):1
  A = Erosion.uplift_map(p, (n, n))
  hm = heatmap!(axis, xs, ys, A; colormap = :grays)
  Colorbar(fig[1, 2], hm)
end

fig, layout, axis = plot()
uplift!(axis, fig, p)
lines!(axis, contour.(0:1/1000:1); color = :red, subdivisions = 1000)
scatter!(axis, contour; color = :green)
lines!(axis, skeleton; color = :blue)
scatter!(axis, skeleton; color = :cyan)
fig

# Tectonic-based erosion

# contour = Patch{BezierCurve,3}(P2[(0.0, 0.0), (0.3, 0.45), (0.3, 0.6), (0.3, 0.7), (0.7, 0.7), (0.8, 0.1), (0.0, 0.0)])
# skeleton = Patch{BezierCurve, 3}(P2[(0.1, 0.1), (0.2, 0.3), (0.3, 0.35), (0.5, 0.4), (0.6, 0.6)])
# data = (; contour, skeleton, uplift = 0.5, radius = 1.5)
# uplift = UpliftPrimitive{Float64}(UPLIFT_PRIMITIVE_ASYMMETRIC_RIDGE, data)
uplift = load_uplift("alpes_noise.png")
model = TectonicBasedErosion(uplift, 100; speed = 1, smooth_factor = 0.00, stream_power = 0.0005, uplift_factor = 0.1, inverse_momentum_power = 4)
# elevation = generate_terrain((256, 256); seed)
elevation = zeros((256, 256))
maps = ErosionMaps(elevation, model)

display(heatmap(maps.uplift))

for i in 1:100
  erode!(maps, model; progress = true)
  display(heatmap(maps.elevation; colormap = :grays))
  # display(heatmap(maps.drainage))
end

# @time erode!(maps, model)
# @profview erode!(maps, model)
