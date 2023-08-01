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
