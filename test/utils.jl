using ProceduralNoise
using Random: default_rng, seed!
using GLMakie

function heatmap(matrix)
  fig, ax, hm = GLMakie.heatmap(matrix; colormap = :inferno)
  Colorbar(fig[:, end+1], hm)
  fig
end

function generate_terrain(resolution; seed = rand(UInt64))
  seed!(default_rng(), seed)
  noise = Fractal{Perlin}((2, 2), octaves = 8)
  noise(resolution)
end

# Requires loading Plots: heatmap; FileIO; ImageIO.
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
