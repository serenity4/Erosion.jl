struct TectonicBasedErosion{D} <: ErosionModel{D}
  duration::Float64
  timestep::Float64
  p::Int
end

TectonicBasedErosion(args...; kwargs...) = TectonicBasedErosion{CPU}(args...; kwargs...)
TectonicBasedErosion{D}(n::Integer; timestep = 0.1, kwargs...) where {D} = TectonicBasedErosion{D}(n * timestep; timestep, kwargs...)
function TectonicBasedErosion{D}(duration::Float64;
    timestep = 0.1,
    p = 4,
  ) where {D}
  TectonicBasedErosion{D}(duration, timestep, p)
end

erode!(elevation, uplift, model::TectonicBasedErosion; progress = false) = erode!(elevation, uplift, model; progress)
execution_state(::TectonicBasedErosion, elevation) = zeros(size(elevation))

function erode!(elevation, drainage_map, model::TectonicBasedErosion, uplift::Union{AbstractMatrix, UpliftPrimitive, UpliftTree}; progress = false)
  uplift = uplift isa Union{UpliftPrimitive, UpliftTree} ? uplift_map(uplift, size(elevation)) : deepcopy(uplift)
  for t in time_range(model)
    progress && print("\r$(Base.text_colors[:green])TectonicBasedErosion$(Base.text_colors[:default]): timestep $(round(t; digits = 1))/$(model.duration)                       ")
    simulate!(elevation, drainage_map, uplift, model)
  end
  ErosionResult(elevation, (; drainage_map, uplift))
end

function simulate!(elevation, drainage_map, uplift, model::TectonicBasedErosion{CPU})
  nx, ny = size(elevation)
  for i in 1:nx
    for j in 1:ny
      simulate!(elevation, drainage_map, uplift, model, GridPoint(i, j), (nx, ny))
    end
  end
end

function simulate!(elevation, drainage_map, uplift, model::TectonicBasedErosion, point, (nx, ny))
  @assert !isnan(drainage_map[point])
  @assert !isnan(elevation[point])
  elevation_change = 0.0
  drainage = 1.0
  for neighbor in neighbors(point, EightNeighbors())
    is_outside_grid(neighbor, (nx, ny)) && continue
    weight = drainage_weight(elevation, point, neighbor, (nx, ny), model.p)
    @assert !isnan(weight)
    @assert 0 ≤ weight ≤ 1
    drainage += weight * drainage_map[neighbor]
  end
  drainage_map[point] = drainage
  elevation[point] += model.timestep * (uplift[point] - drainage + elevation_change)
end

function compute_slope(elevation, point, neighbor)
  Δh = elevation[point] - elevation[neighbor]
  d = norm(neighbor.coords .- point.coords, 2)
  Δh/d
end

function drainage_weight(elevation, point, neighbor, (nx, ny), p)
  denom = 0.0
  for other_neighbor in neighbors(neighbor, EightNeighbors())
    is_outside_grid(other_neighbor, (nx, ny)) && continue
    denom += compute_slope(elevation, neighbor, other_neighbor)^p
  end
  iszero(denom) && return 0.0
  compute_slope(elevation, point, neighbor)^p/denom
end
