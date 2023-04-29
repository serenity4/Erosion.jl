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
  rng::RNG
  seed::UInt64
  terrain_scale::Float64 # in meter per pixel
  evaporation::Float64 # in m²/s
  rain_amount::Float64
  timestep::Float64
end

function HydraulicErosionV2(duration;
    gravity = 9.81,
    rng = default_rng(),
    seed = rand(UInt64),
    terrain_scale = 5.0,
    evaporation = 0.001,
    rain_amount = 0.1,
    timestep = choose_timestep_cfl(terrain_scale),
  )
  HydraulicErosionV2{typeof(rng)}(gravity, duration, rng, seed, terrain_scale, evaporation, rain_amount, timestep)
end

# The timestep is chosen based on the CFL condition 
ESTIMATED_MAX_VELOCITY = 10.0 # in m/s
choose_timestep_cfl(terrain_scale) = (terrain_scale / ESTIMATED_MAX_VELOCITY) ^ 2

time_range(model::HydraulicErosionV2) = zero(model.duration):model.timestep:model.duration

function erode!(terrain, water, water_flow, velocity, sediment, model::HydraulicErosionV2; execution = CPU(model))
  seed!(model.rng, model.seed)
  rainfall = rainfall_map(size(terrain))
  for t in time_range(model)
    add_water!(water, model, rainfall, execution)
    simulate_shallow_water_flow!(water_flow, water, velocity, model, terrain, execution)
    # erode_and_deposit!(terrain, model, velocity, execution)
    # transport_sediment!(sediment, model, velocity, execution)
    # evaporate!(water, model, execution)
  end
  terrain
end

CPU(model::HydraulicErosionV2) = CPU(nothing)
rainfall_map(resolution) = Fractal{Perlin}((2, 2) .^ 4)(resolution)
add_water!(water, model, rainfall, ::CPU) = water .+= rainfall .* model.timestep .* model.rain_amount

function simulate_shallow_water_flow!(water_flow, water, velocity, model, terrain, execution)
  compute_water_flows!(water_flow, model, terrain, water, execution)
  computer_water_height_and_velocity!(water, velocity, model, water_flow, execution)
end

function update_with_loop!(f, A)
  ni, nj = size(A)
  for j in 1:nj, i in 1:ni
    point = GridPoint(i, j)
    A[point] = f(point)
  end
end

function compute_water_flows!(water_flow, model, terrain, water, ::CPU)
  update_with_loop!(point -> water_flows(water_flow, point, model, terrain, water), water_flow)
end

function water_flows(water_flow, point::GridPoint, model::HydraulicErosionV2, terrain, water)
  d = water[point]
  combined_flows = water_flow[point] # every value is a 4-dimensional vector packing flows with left, right, top and bottom cells.
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
    pipe_cross_section = pipe_length * min(d + dadj)
    flow = prev_flow + model.timestep * pipe_cross_section/pipe_length * model.gravity * Δh
    max(zero(flow), flow)
  end
  # Add a scaling factor to avoid having negative water values when taking incoming flows into consideration.
  column_volume = d * column_area(model)
  K = min(1, column_volume / (sum(components) * model.timestep))
  K .* components
end

column_area(model) = model.terrain_scale ^ 2

neighbors(point::GridPoint) = @SVector [point.left, point.right, point.bottom, point.top]
is_outside_grid(point::GridPoint, (ni, nj)) = point[1] in (0, 1 + ni) || point[2] in (0, 1 + nj)

function computer_water_height_and_velocity!(water, velocity, model, water_flow, ::CPU)
  ni, nj = size(water_flow)
  for j in 1:nj, i in 1:ni
    point = GridPoint(i, j)
    prev_height = water[point]
    water[point] = water_height(water, point, model, water_flow)
    velocity[point] = water_velocity(point, model, water, water_flow, prev_height)
  end
end

function water_height(water, point::GridPoint, model::HydraulicErosionV2, water_flow)
  flow_in = sum(ntuple(i -> input_flow(water_flow, point, i), 4))
  flow_out = sum(water_flow[point])
  volume_change = model.timestep * (flow_in - flow_out)
  water[point] + volume_change / column_area(model)
end

function input_flow(water_flow, point, i)
  adj = neighbor(point, i)
  is_outside_grid(adj, size(water_flow)) && return 0.0
  -output_flow(water_flow, adj, reverse_adjacency_index(i))
end

output_flow(water_flow, point, i) = water_flow[point][i]

# Left -> Right, Right -> Left, Bottom -> Top, Top -> Bottom
reverse_adjacency_index(i) = i + 1 - 2iseven(i)

net_flow(water_flow, point, i) = input_flow(water_flow, point, i) - output_flow(water_flow, point, i)
average_flow(water_flow, point, direction) = (net_flow(water_flow, point, direction) + net_flow(water_flow, point, direction + 2)) / 2
average_flow_y(water_flow, point) = input_flow(water_flow, point.left) - output_flow(water_flow, point, 1)

function water_velocity(point::GridPoint, model::HydraulicErosionV2, water, water_flow, prev_height)
  flow_change = ntuple(direction -> average_flow(water_flow, point, direction), 2)
  average_height = (prev_height + water[point] ) / 2
  flow_change ./ (model.terrain_scale * average_height)
end
