using GLMakie
GLMakie.activate!(inline = true)

function plot(; resolution = (1080, 1080), kwargs...)
  fig = Figure(; size = resolution, kwargs...)
  layout = fig[1, 1]
  axis = Axis(layout; aspect = 1)
  fig, layout, axis
end

function plot_heatmap!(fig, axis, matrix; colormap = :inferno)
  hm = heatmap!(axis, matrix; colormap)
  Colorbar(fig[1, 2], hm)
  fig
end

function heatmap(matrix; colormap = :inferno)
  fig, layout, axis = plot()
  plot_heatmap!(fig, axis, matrix; colormap)
end

function erode_and_save(A::Matrix, model::ErosionModel; terrain = "terrain.png", eroded = "eroded.png")
  !isabspath(terrain) && (terrain = joinpath(@__DIR__, terrain))
  !isabspath(eroded) && (eroded = joinpath(@__DIR__, eroded))
  display(heatmap(A))
  save(terrain, A)
  result = @time erode(A, model; progress = true)
  B = result.terrain
  display(heatmap(B))
  low, high = extrema(B)
  if low < 0 || high > 1
    @warn "Terrain values fall outside [0, 1] (min: $low, max: $high), clamping"
    clamp!(B, 0, 1)
  end
  save(eroded, B)
end
