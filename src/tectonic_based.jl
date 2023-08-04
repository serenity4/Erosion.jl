struct TectonicBasedErosion{D,U,P} <: ErosionModel{D}
  uplift::U
  iterations::Integer
  speed::Float64
  stream_power::Float64
  uplift_factor::Float64
  smooth_factor::Float64
  inverse_momentum_power::Val{P}
  min_slope::Float64
end

TectonicBasedErosion(args...; kwargs...) = TectonicBasedErosion{CPU}(args...; kwargs...)
function TectonicBasedErosion{D}(uplift, iterations::Integer = 1;
    speed = 100.0,
    stream_power = 0.0005,
    uplift_factor = 0.01,
    smooth_factor = 0.01,
    inverse_momentum_power = 4,
    min_slope = 0.001,
  ) where {D}
  TectonicBasedErosion{D,typeof(uplift),inverse_momentum_power}(uplift, iterations, speed, stream_power, uplift_factor, smooth_factor, Val(inverse_momentum_power), min_slope)
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
  uplift = model.uplift isa Union{UpliftPrimitive, UpliftTree} ? uplift_map(model.uplift, size(elevation)) : deepcopy(model.uplift)
  ErosionMaps(zeros(size(elevation)), zeros(size(elevation)), elevation, zeros(size(elevation)), uplift, zeros(size(elevation)))
end

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
  Threads.@threads for i in 1:nx
    for j in 1:ny
      simulate!(maps, model, GridPoint(i, j), (nx, ny))
    end
  end
  copyto!(maps.drainage, maps.new_drainage)
  copyto!(maps.elevation, maps.new_elevation)
  copyto!(maps.uplift, maps.new_uplift)
end

function simulate!(maps::ErosionMaps, model::TectonicBasedErosion, point, (nx, ny))
  (; elevation, drainage, uplift) = maps
  @assert !isnan(drainage[point])
  @assert !isnan(elevation[point])
  Δh = laplacian(elevation, point, (nx, ny))
  drainage_value = compute_drainage(drainage, elevation, point, (nx, ny), model)
  maps.new_drainage[point] = drainage_value
  maps.new_elevation[point] += model.speed * (model.uplift_factor * uplift[point] - model.stream_power * drainage_value + model.smooth_factor * Δh)
end

function compute_drainage(drainage, elevation, point, (nx, ny), model)
  value = 0.0
  steepest_slope = 0.0
  steepest_neighbor = point
  for neighbor in neighbors(point, EightNeighbors())
    is_outside_grid(neighbor, (nx, ny)) && continue
    slope = compute_slope(elevation, neighbor, point)
    slope > model.min_slope || continue
    if steepest_slope < slope
      steepest_slope = max(steepest_slope, slope)
      steepest_neighbor = neighbor
    end
    weight = drainage_weight(elevation, point, neighbor, (nx, ny), inverse_momentum_power(model))
    @assert !isnan(weight)
    @assert 0 ≤ weight ≤ 1
    value += weight * drainage[neighbor]
  end
  steepest_slope != 0.0 && (value += norm((1/nx, 1/ny)) * grid_distance(point, steepest_neighbor) * steepest_slope)
  value
end

@inline function compute_slope(elevation, point, neighbor)
  Δh = elevation[point] - elevation[neighbor]
  Δh/grid_distance(point, neighbor)
end

# hypot/norm are slower
# grid_distance(point, neighbor) = sqrt(sum((neighbor.coords .- point.coords) .^ 2))
grid_distance(point, neighbor) = neighbor in neighbors(point, FourNeighbors()) ? 1.0 : sqrt(2)

@inline function drainage_weight(elevation, point, neighbor, (nx, ny), p)
  denom = 0.0
  h = elevation[neighbor]
  # For a given neighbor, how much water flows to `point` depends on how
  # how steep the slopes with its other neighbors are.
  for other_neighbor in neighbors(neighbor, EightNeighbors())
    is_outside_grid(other_neighbor, (nx, ny)) && continue
    elevation[other_neighbor] < h || continue
    denom += compute_slope(elevation, neighbor, other_neighbor)^p
  end
  iszero(denom) && return 0.0
  compute_slope(elevation, neighbor, point)^p/denom
end

function laplacian(elevation, point, (nx, ny))
  i, j = point.coords
  res = zero(eltype(elevation))
  for k in -1:1
    p = GridPoint(i + k, j)
    is_outside_grid(p, (nx, ny)) && continue
    fac = ifelse(iszero(k), 1, -2)
    res += fac * elevation[point]
  end
  for l in -1:1
    p = GridPoint(i, j + l)
    is_outside_grid(p, (nx, ny)) && continue
    fac = ifelse(iszero(l), 1, -2)
    res += fac * elevation[point]
  end
  res
end
