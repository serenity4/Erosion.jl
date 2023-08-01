@enum UpliftOperation begin
  UPLIFT_OPERATION_PRIMITIVE
  UPLIFT_OPERATION_BLEND
  UPLIFT_OPERATION_WARP
  UPLIFT_OPERATION_DIFFERENCE
end

struct UpliftTree
  operation::UpliftOperation
  children::Any
end

@enum UpliftPrimitiveType begin
  UPLIFT_PRIMITIVE_ASYMMETRIC_RIDGE
  UPLIFT_PRIMITIVE_MOUNTAIN_RANGE
  UPLIFT_PRIMITIVE_VALLEY
end

struct UpliftPrimitive{T}
  type::UpliftPrimitiveType
  data::Any
end

Base.broadcastable(x::Union{UpliftTree, UpliftPrimitive}) = Ref(x)

function Base.getproperty(primitive::UpliftPrimitive{T}, name::Symbol) where {T}
  name === :asymmetric_ridge && return primitive.data::@NamedTuple{contour::Patch{BezierCurve,3,Vector{Point{2,T}}}, skeleton::Patch{BezierCurve,3,Vector{Point{2,T}}}, uplift::Float64, radius::Float64}
  name === :mountain_range && return primitive.data::@NamedTuple{radius::Patch{BezierCurve,3}, skeleton::Segment{2}, uplift::Patch{BezierCurve,3}}
  name === :valley && return primitive.data::@NamedTuple{contour::Patch{BezierCurve,3}, skeleton::Patch{BezierCurve,3}}
  getfield(primitive, name)
end

function uplift(primitive::UpliftPrimitive, p)
  @match primitive.type begin
    &UPLIFT_PRIMITIVE_ASYMMETRIC_RIDGE => begin
      (; contour, skeleton, uplift, radius) = primitive.asymmetric_ridge
      p′ = projection(skeleton, p)
      line = Line(p′, p)
      p ≈ p′ && return uplift
      ℐ = filter!(x -> GeometryExperiments.coordinate(line, x) ≥ 0, intersect(line, contour))
      isempty(ℐ) && return 0.0
      m = GeometryExperiments.nearest(PointSet(ℐ), p′)
      d = norm(p′ - p)/norm(p′ - m)
      uplift * uplift_falloff(d, radius)
    end
  end
end

uplift_falloff(distance, radius) = distance > radius ? zero(distance) : (1 - (distance/radius)^2)^3

function uplift_map(primitive, (nx, ny))
  Δx, Δy = 1 ./ ((nx, ny) .- 1)
  grid = [Point2(i, j) for i in 0:Δx:1, j in 0:Δy:1]
  uplift.(primitive, grid)
end
