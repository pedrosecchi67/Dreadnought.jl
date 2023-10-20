cosine_discretization(N::Int) = let θs = collect(LinRange(0.0, π / 2, N))
    @. sin(θs) ^ 2
end

const default_discretization = cosine_discretization(100)

export Selig

include("NACA.jl")

_cubic_spline(x, y) = extrapolate(
    interpolate(x,y,FritschCarlsonMonotonicInterpolation()),
    Flat()
)
