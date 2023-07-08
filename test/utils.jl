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
function erode_and_save(A::Matrix, erosion::HydraulicErosion; terrain = "terrain.png", eroded = "eroded.png")
  !isabspath(terrain) && (terrain = joinpath(@__DIR__, terrain))
  !isabspath(eroded) && (eroded = joinpath(@__DIR__, eroded))
  A = copy(A)
  display(heatmap(A))
  save(terrain, A)
  display(@time erode!(A, erosion; progress = true))
  display(heatmap(A))
  low, high = extrema(A)
  if low < 0 || high > 1
    @warn "Terrain values fall outside [0, 1] (min: $low, max: $high), clamping"
    clamp!(A, 0, 1)
  end
  save(eroded, A)
end
