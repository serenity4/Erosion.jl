struct TectonicBasedErosion{E,U,P} <: ErosionModel{E}
  uplift::U
  iterations::Int
  speed::Float64
  stream_power::Float64
  uplift_factor::Float64
  smooth_factor::Float64
  inverse_momentum_power::Val{P}
  min_slope::Float64
  scale::Tuple{Float64, Float64}
  execution::E
end

TectonicBasedErosion(args...; execution = CPU(), kwargs...) = TectonicBasedErosion{typeof(execution)}(args...; execution, kwargs...)
function TectonicBasedErosion{E}(uplift, iterations::Integer = 1;
    speed = 100.0,
    stream_power = 0.0005,
    uplift_factor = 0.01,
    smooth_factor = 0.01,
    inverse_momentum_power = 4,
    min_slope = 0.001,
    scale = (150_000, 150_000),
    execution,
  ) where {E}
  TectonicBasedErosion{E,typeof(uplift),inverse_momentum_power}(uplift, iterations, speed, stream_power, uplift_factor, smooth_factor, Val(inverse_momentum_power), min_slope, scale, execution)
end

struct ErosionMaps{M}
  drainage::M
  new_drainage::M
  elevation::M
  new_elevation::M
  uplift::M
  new_uplift::M
end

function ErosionMaps(elevation, model::TectonicBasedErosion)
  ErosionMaps{typeof(elevation)}(zeros(size(elevation)), zeros(size(elevation)), elevation .* 1000, zeros(size(elevation)), materialize_uplift(model.uplift, size(elevation)), zeros(size(elevation)))
end

materialize_uplift(uplift::AbstractMatrix, size) = deepcopy(uplift)
materialize_uplift(uplift::Union{UpliftPrimitive, UpliftTree}, size) = uplift_map(uplift, size)
materialize_uplift(uplift::Nothing, size) = zeros(size)

erode(terrain, model::TectonicBasedErosion{CPU}; kwargs...) = erode!(execution_state(model, terrain), model; kwargs...)

execution_state(model::TectonicBasedErosion, elevation) = ErosionMaps(elevation, model)

inverse_momentum_power(::TectonicBasedErosion{<:Any,<:Any,P}) where {P} = P

function erode!(maps::ErosionMaps, model::TectonicBasedErosion; progress = false)
  for i in 1:model.iterations
    progress && print("\r$(Base.text_colors[:green])TectonicBasedErosion$(Base.text_colors[:default]) (speed = $(model.speed)): iteration $i/$(model.iterations)                       ")
    simulate!(maps, model)
  end
  progress && println()
  ErosionResult(maps.elevation, (; maps.drainage, maps.uplift))
end

function simulate!(maps::ErosionMaps, model::TectonicBasedErosion{CPU})
  nx, ny = size(maps.elevation)
  @parallelize model.execution for i in 1:nx
    for j in 1:ny
      simulate!(maps, model, GridPoint(i, j), (nx, ny))
    end
  end
  copyto!(maps.drainage, maps.new_drainage)
  copyto!(maps.elevation, maps.new_elevation)
end

function simulate!(maps::ErosionMaps, model::TectonicBasedErosion, point, (nx, ny))
  (; elevation, drainage, uplift) = maps
  i, j = point.coords

  # Border nodes are fixed to zero (elevation and drainage)
  if i in (1, nx) || j in (1, ny)
    maps.new_elevation[point] = 0.0
    maps.new_drainage[point] = precipitation(model, (nx, ny))
    return
  end

  @assert !isnan(drainage[point])
  @assert !isnan(elevation[point])
  # Δh = laplacian(elevation, point, (nx, ny), model)

  drained = compute_drainage(drainage, elevation, point, (nx, ny), model)
  maps.new_drainage[point] = drained
  downstream = steepest_neighbor_down(elevation, point, (nx, ny), model)
  stream = drained^0.8 * compute_slope(elevation, point, downstream, model, (nx, ny))^2
  height = elevation[point]
  height -= model.speed * model.stream_power * stream
  # height -= model.speed * (model.stream_power * stream - model.smooth_factor * Δh)
  height = max(height, elevation[downstream])
  height += model.speed * model.uplift_factor * uplift[point]
  maps.new_elevation[point] = height
end

precipitation(model, (nx, ny)) = norm(2 .* model.scale ./ (nx, ny))

function compute_drainage(drainage, elevation, point, (nx, ny), model)
  water = precipitation(model, (nx, ny))
  if model.inverse_momentum_power === Val{Inf}()
    water += accumulated_water_steepest(drainage, elevation, point, (nx, ny), model)
  else
    water += accumulated_water(drainage, elevation, point, (nx, ny), model)
  end
  water
end

function accumulated_water(drainage, elevation, point, (nx, ny), model)
  water = 0.0
  for neighbor in neighbors(point, EightNeighbors())
    is_outside_grid(neighbor, (nx, ny)) && continue
    elevation[neighbor] > elevation[point] || continue
    slope = compute_slope(elevation, point, neighbor, model, (nx, ny))
    slope > model.min_slope || continue
    weight = drainage_weight(elevation, point, neighbor, model, (nx, ny), inverse_momentum_power(model))
    @assert !isnan(weight)
    @assert 0 ≤ weight ≤ 1
    water += weight * drainage[neighbor]
  end
  water
end

function accumulated_water_steepest(drainage, elevation, point, (nx, ny), model)
  water = 0.0
  for neighbor in neighbors(point, EightNeighbors())
    is_outside_grid(neighbor, (nx, ny)) && continue
    steepest = steepest_neighbor_down(elevation, neighbor, (nx, ny), model)
    steepest == point && (water += drainage[neighbor])
  end
  water
end

function steepest_neighbor_up(elevation, point, (nx, ny), model)
  candidate = point
  steepest_slope = 0.0
  for neighbor in neighbors(point, EightNeighbors())
    is_outside_grid(neighbor, (nx, ny)) && continue
    elevation[neighbor] > elevation[point] || continue
    slope = compute_slope(elevation, point, neighbor, model, (nx, ny))
    if steepest_slope < slope
      steepest_slope = slope
      candidate = neighbor
    end
  end
  candidate
end
function steepest_neighbor_down(elevation, point, (nx, ny), model)
  candidate = point
  steepest_slope = 0.0
  for neighbor in neighbors(point, EightNeighbors())
    is_outside_grid(neighbor, (nx, ny)) && continue
    elevation[neighbor] < elevation[point] || continue
    slope = compute_slope(elevation, neighbor, point, model, (nx, ny))
    if steepest_slope < slope
      steepest_slope = slope
      candidate = neighbor
    end
  end
  candidate
end

# Positive if point is lower than neighbor (slope vector pointing upwards).
@inline function compute_slope(elevation, from::GridPoint, to::GridPoint, model, (nx, ny))
  from == to && return 0.0
  Δh = elevation[to] - elevation[from]
  Δh/terrain_distance(from, to, model, (nx, ny))
end

terrain_distance(from, to, model, (nx, ny)) = grid_distance(from, to) * hypot((2 .* model.scale ./ (nx, ny))...)
# hypot/norm are slower
grid_distance(from, to) = sqrt(sum((to.coords .- from.coords) .^ 2))

@inline function drainage_weight(elevation, point, neighbor, model, (nx, ny), p)
  denom = 0.0
  h = elevation[neighbor]
  # For a given neighbor, how much water flows to `point` depends on how
  # how steep the slopes with its other neighbors are.
  for other_neighbor in neighbors(neighbor, EightNeighbors())
    is_outside_grid(other_neighbor, (nx, ny)) && continue
    elevation[other_neighbor] < h || continue
    denom += compute_slope(elevation, neighbor, other_neighbor, model, (nx, ny))^p
  end
  iszero(denom) && return 0.0
  compute_slope(elevation, neighbor, point, model, (nx, ny))^p/denom
end

function laplacian(elevation, point, (nx, ny), model)
  i, j = point.coords
  lx = ly = 0.0
  for k in -1:1
    p = GridPoint(i + k, j)
    is_outside_grid(p, (nx, ny)) && continue
    fac = ifelse(iszero(k), -2, 1)
    lx += fac * elevation[p] / (2model.scale[1] / nx)^2
  end
  for l in -1:1
    p = GridPoint(i, j + l)
    is_outside_grid(p, (nx, ny)) && continue
    fac = ifelse(iszero(l), -2, 1)
    ly += fac * elevation[p] / (2model.scale[2] / ny)^2
  end
  lx + ly
end
