abstract type ExecutionType end

struct CPU{T} <: ExecutionType
  data::T
end

struct GPU{T} <: ExecutionType
  data::T
end
