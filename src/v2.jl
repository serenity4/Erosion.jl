"""
A parallelizable physically-based algorithm for hydraulic erosion.

It is based on solving forward in time 2D shallow-water flow equations, adapted to a 2.5D terrain.

Simulation steps are as follows:
1. Water increases due to rainfall or river sources.
2. Flow is simulated with the shallow-water model, from which the velocity field and the water surface may then be computed.
3. Erosion-deposition process is simulated using the velocity field.
4. Suspended sediment is transported by the velocity field.
5. Water decreases due to evaporation.

**Source**: *Fast hydraulic erosion simulation and visualization on GPU* by Mei X., Decaudin P. and Hu B. G, October 2007, 15th Pacific Conference on Computer Graphics and Applications (PG'07) (pp. 47-56), IEEE.
"""
struct HydraulicErosionV2{RNG<:AbstractRNG} <: HydraulicErosion
  gravity::Float64 # in m/s²
  duration::Int64 # in seconds
  dissolution_constant::Float64
  deposition_constant::Float64
  rng::RNG
  seed::UInt64
  terrain_scale::Float64 # in meters per pixel
  height_scale::Float64 # in meters per unit
  evaporation::Float64 # in m²/s
  rain_amount::Float64
  sediment_transport_capacity_factor::Float64
  minimum_sediment_transport_capacity::Float64
  timestep::Float64
end

function HydraulicErosionV2(duration;
    gravity = 9.81,
    dissolution_constant = 1.0,
    deposition_constant = 1.0,
    rng = default_rng(),
    seed = rand(UInt64),
    terrain_scale = 5.0,
    evaporation = 0.1,
    height_scale = 700.0,
    rain_amount = 0.02,
    sediment_transport_capacity_factor = 1.0,
    minimum_sediment_transport_capacity = 0.0001,
    timestep = choose_timestep_cfl(terrain_scale),
  )
  HydraulicErosionV2{typeof(rng)}(gravity, duration, dissolution_constant, deposition_constant, rng, seed, terrain_scale, height_scale, evaporation, rain_amount, sediment_transport_capacity_factor, minimum_sediment_transport_capacity, timestep)
end

# The timestep is chosen based on the CFL condition 
ESTIMATED_MAX_VELOCITY = 30.0 # in m/s
choose_timestep_cfl(terrain_scale) = terrain_scale / ESTIMATED_MAX_VELOCITY

time_range(model::HydraulicErosionV2) = zero(model.duration):model.timestep:model.duration

erode!(terrain, model::HydraulicErosionV2; execution = CPU(model), progress = false) = erode!(terrain, initialize_maps(typeof(model), terrain)..., model; execution, progress)

function initialize_maps(::Type{<:HydraulicErosionV2}, terrain)
  water = zeros(size(terrain))
  water_flow = zeros(SVector{4,Float64}, size(terrain))
  velocity = zeros(SVector{2,Float64}, size(terrain))
  sediment = zeros(size(terrain))
  water, water_flow, velocity, sediment
end

function erode!(terrain, water, water_flow, velocity, sediment, model::HydraulicErosionV2; execution = CPU(model), progress = false)
  seed!(model.rng, model.seed)
  rainfall = rainfall_map(size(terrain))
  for t in time_range(model)
    progress && print("\r$(Base.text_colors[:green])HydraulicErosionV2$(Base.text_colors[:default]): timestep $(round(t; digits = 1))/$(model.duration)                       ")
    erode!(terrain, water, water_flow, velocity, sediment, model, t, rainfall; execution)
  end
  terrain
end

mean(A) = sum(x -> x^2, A)/length(A)
stats(A) = string("minimum = ", minimum(A), ", mean = ", mean(A), ", maximum = ", maximum(A))

function erode!(terrain, water, water_flow, velocity, sediment, model::HydraulicErosionV2, t::Number, rainfall = rainfall_map(size(terrain)); execution = CPU(model))
  add_water!(water, model, rainfall, execution)
  # @show stats(water)
  simulate_shallow_water_flow!(water_flow, water, velocity, model, terrain, execution)
  # @show stats(norm.(water_flow)) stats(norm.(velocity)) stats(water)
  erode_and_deposit!(terrain, sediment, model, velocity, execution)
  transport_sediment!(sediment, model, velocity, execution)
  evaporate!(water, model, execution)
  # @show stats(terrain) stats(sediment)
end

CPU(model::HydraulicErosionV2) = CPU(nothing)
rainfall_map(resolution) = Fractal{Perlin}((2, 2) .^ 4)(resolution)
add_water!(water, model, rainfall, ::CPU) = water .+= rainfall .* model.timestep .* model.rain_amount
evaporate!(water, model, ::CPU) = water .*= (1 - model.evaporation * model.timestep)

function simulate_shallow_water_flow!(water_flow, water, velocity, model, terrain, execution)
  compute_water_flows!(water_flow, model, terrain, water, execution)
  computer_water_height_and_velocity!(water, velocity, model, water_flow, execution)
end

function update_with_loop!(f, A)
  ni, nj = size(A)
  # Threads.@threads for j in 1:nj
  for j in 1:nj
    for i in 1:ni
      point = GridPoint(i, j)
      A[point] = f(point)
    end
  end
end

function compute_water_flows!(water_flow, model, terrain, water, ::CPU)
  update_with_loop!(point -> water_flows(water_flow, point, model, terrain, water), water_flow)
end

function water_flows(water_flow, point::GridPoint, model::HydraulicErosionV2, terrain, water)
  d = water[point]
  # Consider there to be no water at all (and subsequently no flow) when water height is < 0.01 mm
  combined_flows = water_flow[point] # every value is a 4-dimensional vector packing flows with left, right, bottom and top cells.
  isapprox(d, zero(d); atol = 1e-5) && return @SVector zeros(Float64, 4)
  components = map(enumerate(neighbors(point))) do (i, adj)
    is_outside_grid(adj, size(water_flow)) && return 0.0
    prev_flow = combined_flows[i]
    dadj = water[adj]
    Δh = terrain[point] + d - (terrain[adj] + dadj)
    # Let's take half of the distance between two pixels.
    pipe_length = 0.5 * model.terrain_scale
    # Let's take that same distance for the pipe width, and multiply by the smallest water height.
    # Yes, that means pipes are variable depending on the water level.
    # We'd want large pipes for an ocean, and tiny pipes for shallow streams.
    pipe_cross_section = pipe_length * min(d, dadj)
    flow = prev_flow + model.timestep * pipe_cross_section/pipe_length * model.gravity * Δh
    max(zero(flow), flow)
  end
  all(isapprox(x, zero(x), atol = 1e-5) for x in components) && return @SVector zeros(Float64, 4)
  # Add a scaling factor to avoid having negative water values when taking incoming flows into consideration.
  K = min(1, d * model.terrain_scale^2 / (sum(components) * model.timestep))
  K .* components
end

neighbors(point::GridPoint) = @SVector [point.left, point.right, point.bottom, point.top]

function computer_water_height_and_velocity!(water, velocity, model, water_flow, ::CPU)
  ni, nj = size(water_flow)
  # Threads.@threads for j in 1:nj
  for j in 1:nj
    for i in 1:ni
      point = GridPoint(i, j)
      prev_height = water[point]
      water[point] = water_height(water, point, model, water_flow)
      velocity[point] = water_velocity(point, model, water, water_flow, prev_height)
    end
  end
end

function water_height(water, point::GridPoint, model::HydraulicErosionV2, water_flow)
  flow_in = sum(ntuple(i -> input_flow(water_flow, point, i), 4))
  flow_out = sum(water_flow[point])
  @assert all(≥(0), flow_in)
  @assert all(≥(0), flow_out)
  volume_change = model.timestep * (flow_in - flow_out)
  water[point] + volume_change / model.terrain_scale^2
end

function input_flow(water_flow, point, i)
  adj = neighbor(point, i)
  is_outside_grid(adj, size(water_flow)) && return 0.0
  output_flow(water_flow, adj, reverse_adjacency_index(i))
end

output_flow(water_flow, point, i) = water_flow[point][i]

# Left -> Right, Right -> Left, Bottom -> Top, Top -> Bottom
reverse_adjacency_index(i) = i + 1 - 2iseven(i)

net_flow(water_flow, point, i) = input_flow(water_flow, point, i) - output_flow(water_flow, point, i)
average_flow(water_flow, point, direction) = (net_flow(water_flow, point, direction) + net_flow(water_flow, point, direction + 2)) / 2
average_flow_y(water_flow, point) = input_flow(water_flow, point.left) - output_flow(water_flow, point, 1)

function water_velocity(point::GridPoint, model::HydraulicErosionV2, water, water_flow, prev_height)
  average_height = (prev_height + water[point]) / 2
  isapprox(average_height, zero(average_height); atol = 1e-5) && return @SVector zeros(Float64, 2)
  flow_change = ntuple(direction -> average_flow(water_flow, point, direction), 2)
  flow_change ./ (model.terrain_scale * average_height)
end

function erode_and_deposit!(terrain, sediment, model, velocity, ::CPU)
  ni, nj = size(terrain)
  # Threads.@threads for j in 1:nj
  for j in 1:nj
    for i in 1:ni
      erode_and_deposit!(terrain, sediment, GridPoint(i, j), model, velocity)
    end
  end
end

function erode_and_deposit!(terrain, sediment, point::GridPoint, model::HydraulicErosionV2, velocity)
  α = tilt_angle(terrain, point, model)
  sediment_transport_capacity = norm(velocity[point]) * model.sediment_transport_capacity_factor * sin(α)
  sediment_transport_capacity = max(sediment_transport_capacity, model.minimum_sediment_transport_capacity)
  if sediment_transport_capacity > sediment[point]
    dissolved_sediment = model.dissolution_constant * (sediment_transport_capacity - sediment[point])
    terrain[point] -= dissolved_sediment
    sediment[point] += dissolved_sediment
  else
    deposited_sediment = model.deposition_constant * (sediment[point] - sediment_transport_capacity)
    terrain[point] += deposited_sediment
    sediment[point] -= deposited_sediment
  end
end

function tilt_angle(terrain, point, model)
  gradient = estimate_gradient(terrain, point, size(terrain))
  Δh = norm(gradient)
  atan(Δh, sign(Δh))
  # atan(Δh * model.height_scale / model.terrain_scale ^ 2, sign(Δh))
end

transport_sediment!(sediment, model, velocity, ::CPU) = update_with_loop!(point -> transport_sediment(sediment, point, model, velocity), sediment)

function transport_sediment(sediment, point::GridPoint, model::HydraulicErosionV2, velocity)
  prev = point[] .+ velocity[point] .* model.timestep
  bounds = size(sediment)
  is_outside_grid(prev, bounds) && return sediment[point]
  interpolate_bilinear(sediment, prev, Cell(prev, bounds))
end
