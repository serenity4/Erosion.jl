using Random: default_rng, seed!

function generate_terrain(resolution; seed = rand(UInt64))
  scale = (2, 2) .^ 4
  coords = Tuple(0.0:1/r:(1.0 - 1.0/r) for r in resolution)
  grid = collect(Iterators.product(coords...))
  seed!(default_rng(), seed)
  perlin = Perlin(scale)
  noise = Fractal(perlin, octaves = 8)
  terrain = noise.(grid)
  remap.(terrain, Ref((minimum(terrain), maximum(terrain))), Ref((0.0, 1.0)))
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
