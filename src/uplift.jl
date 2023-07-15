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
      _, p′ = project(skeleton, p)
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
