module Erosion

using Random: seed!, AbstractRNG, default_rng

Base.@kwdef struct HydraulicErosion{RNG<:AbstractRNG}
  gravity::Float64 = 9.81
  iterations::Int64 = 10000
  rng::RNG = default_rng()
  seed::Int64 = 0
  timestep::Float64 = 1.0
  # Minimum slope so that even flat terrains get eroded.
  min_slope::Float64 = 0.01
  erosion_factor::Float64 = 1.0
  deposition_factor::Float64 = 1.0
  evaporation::Float64 = 0.05
  terrain_size::Tuple{Float64,Float64} = (1.0, 1.0)
  "Inertia in terms of direction taken by the droplet. Speed is not affected."
  droplet_inertia::Float64 = 0.5
  "Minimum amount of sediment that is deposited (when possible) at every step."
  droplet_min_deposited::Float64 = 0.01
  droplet_capacity::Float64 = 1.0
  droplet_effect_radius::Float64 = 0.01
  droplet_max_steps::Int = 1000
  droplet_min_speed::Float64 = 0.05
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

nearest((x, y)) = Int.(round.((x, y)))

function Cell((x, y))
  (i, j) = Int.(floor.((x, y)))
  Cell((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1))
end

DropletSampler(terrain_size, i) = DropletSampler(terrain_size)

function Base.rand(rng::AbstractRNG, sampler::DropletSampler)
  ni, nj = sampler.terrain_size
  position = (rand(rng, 1:(ni - 1)) + rand(rng), rand(rng, 1:(nj - 1)) + rand(rng))
  Droplet(position, (0.0, 0.0), 0.0, 0.0, 1.0)
end

function erode!(terrain, model::HydraulicErosion)
  seed!(model.rng, model.seed)
  codes = DropletResult[]
  for i in 1:model.iterations
    droplet = rand(model.rng, DropletSampler(size(terrain), i))
    code = simulate!(terrain, model, droplet)
    !in(code, codes) && push!(codes, code)
  end
  codes
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
  BASSIN = 2
  ESCAPED = 3
end

"Went over the iteration limit defined `droplet_max_steps`."
REACHED_ITERATION_LIMIT
"Evaporated all its water."
EVAPORATED
"Reached a pit that the droplet could not completely fill with sediment."
BASSIN

function simulate!(terrain, model::HydraulicErosion, droplet::Droplet)
  units_per_pixel = model.terrain_size ./ size(terrain)
  effect_radius = model.droplet_effect_radius ./ units_per_pixel
  droplet_nx, droplet_ny = neighborhood(effect_radius)
  falloff_weights = Float64[elliptic_falloff((0.0, 0.0), (i, j), effect_radius) for i in droplet_nx, j in droplet_ny]
  falloff_weights ./= sum(falloff_weights)
  for _ in 1:model.droplet_max_steps
    # Evaporation.
    water_amount = droplet.water_amount * (1 - model.evaporation)
    isapprox(water_amount, zero(water_amount); atol = 1e-4) && return EVAPORATED

    # Compute new position and direction.
    cell = Cell(droplet.position)
    h₀ = interpolate_bilinear(terrain, droplet.position, cell)
    slope = estimate_gradient(terrain, droplet.position, cell)
    # Compute new droplet data.
    direction = droplet.direction .* model.droplet_inertia .- slope .* (1 - model.droplet_inertia)
    iszero(norm(direction)) && (direction = (rand(), rand()))
    direction = normalize(direction)
    # The position update is designed to be spatially constant, sacrificing
    # temporal consistency with respect to the speed of the particle.
    position = droplet.position .+ direction
    !in_terrain(terrain, position) && return ESCAPED
    cell = Cell(position)
    h = interpolate_bilinear(terrain, position, cell)
    Δh = h - h₀

    # Optionally fill bassin if we reached a dip.
    if Δh > 0
      # We went up.
      # We don't go up (obviously), and instead deposit sediment in the bassin.
      # If there is enough sediment to fill the bassin, then we keep going downward
      # with the newly computed height.
      # TODO
      # Δh = ...
    end
    if Δh > 0
      # Couldn't fill the bassin.
      return BASSIN
    end

    # We are going down.
    α = norm(units_per_pixel)/-Δh # inclination, 0 = flat terrain, 1 = straight cliff)
    # @assert 0 ≤ α ≤ 1 "Expected inclination to be between [0, 1], got α = $α"
    # We make sure we always erode even when the terrain is flat (`max(...)`).
    # This could be changed to never erode when the terrain is flat.
    speed = max(α * model.gravity, model.droplet_min_speed)

    # Sedimentation.
    ## Deposition.
    ## Some is deposited by evaporation, some is diffused along the way.
    sediment_capacity = water_amount * model.droplet_capacity
    (; carried_sediment) = droplet
    evaporated_sediment = min(0, carried_sediment - sediment_capacity)
    carried_sediment -= evaporated_sediment
    diffused_sediment = min(max(model.droplet_min_deposited, (1 - α)^2) / speed * model.erosion_factor * carried_sediment, carried_sediment)
    carried_sediment -= diffused_sediment
    deposited_sediment = evaporated_sediment + diffused_sediment

    ## Ablation by erosion.
    eroded_sediment = (sediment_capacity - carried_sediment) * min(speed * model.erosion_factor, 1)
    carried_sediment += eroded_sediment
    carried_sediment = min(carried_sediment, sediment_capacity)

    ## Update neighboring terrain.
    from_pixel = nearest(droplet.position)
    to_pixel = nearest(position)
    for i in droplet_nx, j in droplet_ny
      u, v = (i, j) .+ from_pixel
      in_terrain(terrain, (u, v)) && (terrain[u, v] -= falloff_weights[1 + i - droplet_nx.start, 1 + j - droplet_ny.start] * eroded_sediment)
      u, v = (i, j) .+ to_pixel
      in_terrain(terrain, (u, v)) && (terrain[u, v] += falloff_weights[1 + i - droplet_nx.start, 1 + j - droplet_ny.start] * deposited_sediment)
    end

    droplet = Droplet(position, direction, speed, carried_sediment, water_amount)
  end
  REACHED_ITERATION_LIMIT
end

in_terrain(terrain, coords) = all(1 ≤ c ≤ n for (c, n) in zip(coords, size(terrain)))

# Later we can compute this factor in a more dynamic way to control where erosion should occur.
erosion_factor(model::HydraulicErosion, (x, y, z)) = model.erosion_factor

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
