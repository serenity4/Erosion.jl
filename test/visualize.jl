using ImageTransformations

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
nx, ny = (256, 256)
uplift = load_uplift("alpes_noise.png")
model = TectonicBasedErosion(uplift, 1; speed = 100, smooth_factor = 0.00, stream_power = 0.0005, uplift_factor = 0.1, inverse_momentum_power = 4)
# elevation = generate_terrain((256, 256); seed)
elevation = imresize(load_elevation("fractal.png"), (nx, ny)) .* 10
maps = ErosionMaps(elevation, model)

display(heatmap(maps.uplift))
display(heatmap(elevation))

for i in 1:100
  erode!(maps, model; progress = true)
  display(heatmap(maps.elevation; colormap = :grays))
  # display(heatmap(maps.drainage))
end

# @time erode!(maps, model)
# @profview erode!(maps, model)

# Drainage networks

nx, ny = (256, 256)
elevation = generate_terrain((nx, ny); seed)
elevation = [(2x - 1)^2 + (2y - 1)^2 for x in 0:1/(nx - 1):1, y in 0:1/(ny - 1):1]
elevation = ones(nx, ny)
elevation[nx ÷ 2, ny ÷ 2] = 0.0
elevation = load_uplift("alpes_noise.png")
elevation = imresize(load_elevation("fractal.png"), (nx, ny))
elevation = imresize(load_elevation("mountains.png"), (nx, ny))
elevation = imresize(load_elevation("grand_mountain.png"), (nx, ny))
display(heatmap(elevation))

model = TectonicBasedErosion(nothing; inverse_momentum_power = 4)
drainage = zeros((nx, ny))
new_drainage = zeros((nx, ny))
for i in 1:100
  Threads.@threads for i in 1:nx
    for j in 1:ny
      point = GridPoint(i, j)
      new_drainage[point] = compute_drainage(drainage, elevation, point, (nx, ny), model)
    end
  end
  copyto!(drainage, new_drainage)
  i % 10 == 1 && display(heatmap(drainage; colormap = :grays))
  # display(heatmap(drainage; colormap = :grays))
end
display(heatmap(drainage; colormap = :grays))

fig, layout, axis = plot(; resolution = (1920, 1080))
axis2 = Axis(fig[1, 2]; aspect = 1)
heatmap!(axis, elevation; colormap = :grays)
heatmap!(axis2, drainage; colormap = :grays)
fig

using Colors, ImageShow
using SPIRV.MathFunctions

A = convert(Matrix{RGB{Float64}}, remap.(elevation, extrema(elevation)..., 0.0, 1.0))
alpha = remap.(drainage, extrema(drainage)..., 0.0, 1.0)
blend(src, dst, alpha) = src .* (1 .- alpha) .+ dst .* alpha
B = blend(A, RGB{Float64}(1, 0, 0), alpha)
save("elevation_with_drainage.png", B)
save("elevation_with_drainage_2.png", B)
