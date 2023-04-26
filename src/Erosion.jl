module Erosion

using Random: seed!, AbstractRNG, default_rng

struct HydraulicErosion{RNG<:AbstractRNG}
  gravity::Float64
  iterations::Int64
  rng::RNG
  seed::UInt64
  # Minimum slope so that even flat terrains get eroded.
  min_slope::Float64
  # Later we can compute this factor in a more dynamic way to control where erosion should occur.
  erosion_factor::Float64
  deposition_factor::Float64
  evaporation::Float64
  terrain_size::Tuple{Float64,Float64}
  "Inertia in terms of direction taken by the droplet. Speed is not affected."
  droplet_inertia::Float64
  "Minimum amount of sediment that is deposited (when possible) at every step."
  droplet_min_deposited::Float64
  droplet_capacity::Float64
  droplet_effect_radius::Float64
  droplet_max_steps::Int
  droplet_min_speed::Float64
end

HydraulicErosion(gravity, iterations, rng, args...) = HydraulicErosion{typeof(rng)}(gravity, iterations, rng, args...)

function HydraulicErosion(;
    gravity = 9.81,
    iterations = 10000,
    rng = default_rng(),
    seed = rand(UInt64),
    min_slope = 0.01,
    erosion_factor = 1.0,
    deposition_factor = 1.0,
    evaporation = 0.05,
    terrain_size = (1.0, 1.0),
    droplet_inertia = 0.4,
    droplet_min_deposited = 0.01,
    droplet_capacity = 1.0,
    droplet_effect_radius = 0.01,
    droplet_max_steps = 1000,
    droplet_min_speed = 0.05,
  )
  HydraulicErosion(gravity, iterations, rng, seed, min_slope, erosion_factor, deposition_factor, evaporation, terrain_size, droplet_inertia, droplet_min_deposited, droplet_capacity, droplet_effect_radius, droplet_max_steps, droplet_min_speed)
end

struct Droplet
  position::Tuple{Float64,Float64}
  direction::Tuple{Float64,Float64}
  speed::Float64
  carried_sediment::Float64
  water_amount::Float64
end

struct DropletSampler
  terrain_size::Tuple{Int64,Int64}
end

DropletSampler(terrain_size, i) = DropletSampler(terrain_size)

function Base.rand(rng::AbstractRNG, sampler::DropletSampler)
  ni, nj = sampler.terrain_size
  position = (rand(rng, 1:(ni - 1)) + rand(rng), rand(rng, 1:(nj - 1)) + rand(rng))
  Droplet(position, (0.0, 0.0), 0.0, 0.0, 1.0)
end

struct GridPosition
  coords::Tuple{Int,Int}
end

Base.iterate(pos::GridPosition, args...) = iterate(pos.coords, args...)
Base.convert(::Type{GridPosition}, coords::Tuple{Int,Int}) = GridPosition(coords)

struct Cell
  bottom_left::GridPosition
  bottom_right::GridPosition
  top_right::GridPosition
  top_left::GridPosition
end

Base.getindex(A::AbstractArray, gpos::GridPosition) = A[gpos.coords...]
Base.setindex!(A::AbstractArray, value, gpos::GridPosition) = A[gpos.coords...] = value

nearest((x, y)) = Int.(round.((x, y)))

function Cell((x, y))
  (i, j) = Int.(floor.((x, y)))
  Cell((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1))
end

function bilinear_weights(cell::Cell, (x, y))
  (cx, cy) = cell.bottom_left
  w1 = (-(x - (1 + cx)) * -(y - (1 + cy)))
  w2 = ((x - cx) * -(y - (1 + cy)))
  w3 = (-(x - (1 + cx)) * (y - cy))
  w4 = ((x - cx) * (y - cy))
  (w1, w2, w3, w4)
end

corners(cell::Cell) = (cell.bottom_left, cell.bottom_right, cell.top_left, cell.top_right)

struct ErosionMetrics
  reached_iteration_limit::Float64
  evaporated::Float64
  basin::Float64
  escaped::Float64
end

function ErosionMetrics(codes)
  n = length(codes)
  results = (REACHED_ITERATION_LIMIT, EVAPORATED, BASIN, ESCAPED)
  ErosionMetrics(ntuple(i -> count(==(results[i]), codes)/n, length(results))...)
end

function Base.show(io::IO, metrics::ErosionMetrics)
  print(io, ErosionMetrics, '(')
  for (i, name) in enumerate(fieldnames(ErosionMetrics))
    i > 1 && print(io, ", ")
    val = getproperty(metrics, name)
    print(io, name, " = ", round(100val; digits=2), '%')
  end
  print(io, ')')
end

function erode!(terrain, model::HydraulicErosion)
  seed!(model.rng, model.seed)
  codes = DropletResult[]
  units_per_pixel = model.terrain_size ./ size(terrain)
  effect_radius = model.droplet_effect_radius ./ units_per_pixel
  window = neighborhood(effect_radius)
  falloff_weights = Float64[elliptic_falloff((0.0, 0.0), (i, j), effect_radius) for i in window[1], j in window[2]]
  falloff_weights ./= sum(falloff_weights)
  for i in 1:model.iterations
    droplet = rand(model.rng, DropletSampler(size(terrain), i))
    code = simulate!(terrain, model, droplet, units_per_pixel, effect_radius, window, falloff_weights)
    push!(codes, code)
  end
  ErosionMetrics(codes)
end

function interpolate_bilinear(A, (x, y), cell::Cell)
  cx, cy = cell.bottom_left
  nx0 = lerp(A[cell.bottom_left], A[cell.bottom_right], x - cx)
  nx1 = lerp(A[cell.top_left], A[cell.top_right], x - cx)
  lerp(nx0, nx1, y - cy)
end

lerp(x, y, w) = x * (1 - w) + y * w

function estimate_gradient(A, (x, y), cell::Cell)
  cx, cy = cell.bottom_left
  gx = lerp(A[cell.bottom_right] - A[cell.bottom_left], A[cell.top_right] - A[cell.top_left], x - cx)
  gy = lerp(A[cell.top_right] - A[cell.bottom_right], A[cell.top_left] - A[cell.bottom_left], y - cy)
  (gx, gy)
end

norm((x, y)) = hypot(x, y)
normalize(direction) = direction ./ norm(direction)

@enum DropletResult begin
  REACHED_ITERATION_LIMIT = -1
  EVAPORATED = 1
  BASIN = 2
  ESCAPED = 3
end

"Went over the iteration limit defined `droplet_max_steps`."
REACHED_ITERATION_LIMIT
"Evaporated all its water."
EVAPORATED
"Reached a pit that the droplet could not completely fill with sediment."
BASIN

function simulate!(terrain, model::HydraulicErosion, droplet::Droplet, units_per_pixel, effect_radius, (droplet_nx, droplet_ny), falloff_weights)
  for _ in 1:model.droplet_max_steps
    # Evaporation.
    water_amount = droplet.water_amount * (1 - model.evaporation)
    isapprox(water_amount, zero(water_amount); atol = 1e-4) && return EVAPORATED

    # Compute new position and direction.
    from_cell = Cell(droplet.position)
    h₀ = interpolate_bilinear(terrain, droplet.position, from_cell)
    slope = estimate_gradient(terrain, droplet.position, from_cell)
    # Compute new droplet data.
    direction = droplet.direction .* model.droplet_inertia .- slope .* (1 - model.droplet_inertia)
    iszero(norm(direction)) && (direction = (rand(), rand()))
    direction = normalize(direction)
    # The position update is designed to be spatially constant, sacrificing
    # temporal consistency with respect to the speed of the particle.
    position = droplet.position .+ direction
    !in_terrain(terrain, position) && return ESCAPED
    to_cell = Cell(position)
    h = interpolate_bilinear(terrain, position, to_cell)
    Δh = h - h₀

    # Optionally fill basin if we reached a dip.
    if Δh > 0
      # We went up.
      # We don't go up (obviously), and instead deposit sediment in the basin.
      # If there is enough sediment to fill the basin, then we keep iterating.
      # Don't fill more than the height difference, otherwise we would be able to dig holes.
      deposited_sediment = min(Δh, droplet.carried_sediment)
      fill_with_sediment!(terrain, droplet.position, Cell(droplet.position), deposited_sediment)
      # Couldn't fill the basin.
      deposited_sediment == droplet.carried_sediment && return BASIN
      droplet = Droplet(droplet.position, droplet.direction, droplet.speed, 0.0, droplet.water_amount)
      continue
    end

    # We are going down.
    # Inclination, 0 = flat terrain, 1 = straight cliff.
    α = cos(atan(-Δh, norm(units_per_pixel))) # we move 1 pixel-length at a time.
    # We make sure we always erode even when the terrain is flat (`max(...)`).
    # This could be changed to never erode when the terrain is flat.
    speed = max(α * model.gravity, model.droplet_min_speed)

    # Sedimentation.
    ## Deposition.
    ## Some is deposited by evaporation, some is diffused along the way.
    sediment_capacity = water_amount * model.droplet_capacity
    (; carried_sediment) = droplet
    evaporated_sediment = max(0, carried_sediment - sediment_capacity)
    carried_sediment -= evaporated_sediment
    diffused_sediment = min(max(model.droplet_min_deposited, (1 - α)^2) / speed * model.deposition_factor * carried_sediment, carried_sediment)
    carried_sediment -= diffused_sediment
    deposited_sediment = min(evaporated_sediment + diffused_sediment, -Δh)

    ## Ablation by erosion.
    erosion_strength = (droplet.water_amount * model.droplet_capacity - droplet.carried_sediment) * droplet.speed * model.erosion_factor

    @assert deposited_sediment ≥ 0
    @assert carried_sediment ≥ 0

    ## Update neighboring terrain.
    fill_with_sediment!(terrain, position, to_cell, deposited_sediment)
    from_pixel = nearest(droplet.position)

    eroded_sediment = 0.0
    for i in droplet_nx, j in droplet_ny
      u, v = (i, j) .+ from_pixel
      in_terrain(terrain, (u, v)) || continue
      weight = falloff_weights[1 + i - droplet_nx.start, 1 + j - droplet_ny.start]
      hloc = terrain[u, v]
      # Don't erode areas that are at a lower height than the current point.
      h < hloc || continue
      Δhloc = h - hloc
      eroded = erosion_strength * weight * sqrt(-Δhloc)
      eroded = min(eroded, -Δhloc)
      terrain[u, v] -= eroded
      @assert terrain[u, v] ≥ 0
      eroded_sediment += eroded
    end

    @assert eroded_sediment ≥ 0
    carried_sediment += eroded_sediment
    carried_sediment = min(carried_sediment, sediment_capacity)

    droplet = Droplet(position, direction, speed, carried_sediment, water_amount)
  end
  REACHED_ITERATION_LIMIT
end

function fill_with_sediment!(terrain, position, cell::Cell, deposited_sediment)
  for (corner, weight) in zip(corners(cell), bilinear_weights(cell, position))
    fill_with_sediment!(terrain, corner, deposited_sediment * weight)
  end
end

function fill_with_sediment!(terrain, pixel, amount)
  @assert amount ≥ 0
  terrain[pixel] += amount
end

in_terrain(terrain, coords) = all(1 ≤ c ≤ n for (c, n) in zip(coords, size(terrain)))

function neighborhood(radius)
  upper_bound = Int.(ceil.(radius))
  (:).(.-(upper_bound), 1, upper_bound)
end

elliptic_distance(p1, p2, radius) = norm((p2 .- p1) ./ radius)

function elliptic_falloff(origin, position, radius)
  d = elliptic_distance(origin, position, radius)
  max(0.0, 1.0 - d)
end

export HydraulicErosion, erode!

end
