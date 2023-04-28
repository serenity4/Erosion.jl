using Random: default_rng, seed!

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
  save(eroded, A)
end
