abstract type ExecutionType end

struct CPU{T} <: ExecutionType
  data::T
end

struct GPU{T} <: ExecutionType
  data::T
end

abstract type ErosionModel{ExecutionType} end

"""
    erode(terrain, model::ErosionModel)

Erode the provided terrain, returning an [`ErosionResult`](@ref).

The data contained in the `ErosionResult` depends on the model.
"""
function erode end

"""
Result holding an eroded `terrain` and model-specific `data` which may be
of further interest e.g. for rendering or running additional algorithms.
"""
struct ErosionResult{T,D}
  terrain::T
  data::D
end

"""
    execution_state(model::ErosionModel, terrain)

State required for the model to execute, involving e.g. auxiliary maps.
"""
function execution_state end

erode(terrain, model::ErosionModel{CPU}; kwargs...) = erode!(copy(terrain), execution_state(model, terrain), model; kwargs...)
