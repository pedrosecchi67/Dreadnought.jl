export Airfoil

include("airfoil_management.jl")

"""
```
    struct Airfoil
        thickness::AbstractInterpolation
        camber::AbstractInterpolation
        LE_flap_chord::Real
        LE_flap_deflection::Real
        TE_flap_chord::Real
        TE_flap_deflection::Real
        forced_separation::Real
    end
```

Struct to define an airfoil
"""
struct Airfoil
    thickness::AbstractInterpolation
    camber::AbstractInterpolation
    LE_flap_chord::Real
    LE_flap_deflection::Real
    TE_flap_chord::Real
    TE_flap_deflection::Real
    forced_separation::Real
end

"""
```
    Airfoil(
        xs::AbstractVector, ts::AbstractVector, ys::AbstractVector;
        LE_flap_chord::Real = 0.0,
        LE_flap_deflection::Real = 0.0,
        TE_flap_chord::Real = 0.0,
        TE_flap_deflection::Real = 0.0,
        forced_separation::Real = 1.0,
    ) = Airfoil(
        _cubic_spline(xs, ts),
        _cubic_spline(xs, ys),
        LE_flap_chord, LE_flap_deflection,
        TE_flap_chord, TE_flap_deflection,
        forced_separation
    )
```

Define an airfoil from vectors describing 
x axis control points, local thicknesses and local camber, 
respectively
"""
Airfoil(
    xs::AbstractVector, ts::AbstractVector, ys::AbstractVector;
    LE_flap_chord::Real = 0.0,
    LE_flap_deflection::Real = 0.0,
    TE_flap_chord::Real = 0.0,
    TE_flap_deflection::Real = 0.0,
    forced_separation::Real = 1.0,
) = Airfoil(
    _cubic_spline(xs, ts),
    _cubic_spline(xs, ys),
    LE_flap_chord, LE_flap_deflection,
    TE_flap_chord, TE_flap_deflection,
    forced_separation
)

"""
```
    function Airfoil(
        pts::AbstractMatrix; 
        N::Int = 40,
        LE_flap_chord::Real = 0.0,
        LE_flap_deflection::Real = 0.0,
        TE_flap_chord::Real = 0.0,
        TE_flap_deflection::Real = 0.0,
        forced_separation::Real = 1.0,
    )
```

Build airfoil from matrix of Selig-format data points
"""
function Airfoil(
    pts::AbstractMatrix; 
    N::Int = 40,
    LE_flap_chord::Real = 0.0,
    LE_flap_deflection::Real = 0.0,
    TE_flap_chord::Real = 0.0,
    TE_flap_deflection::Real = 0.0,
    forced_separation::Real = 1.0,
)

    x = @view pts[:, 1]
    y = @view pts[:, 2]

    ile = argmin(x)

    xup = x[ile:-1:1]
    yup = y[ile:-1:1]

    xlow = x[ile:end]
    ylow = y[ile:end]

    lower, upper = (
        _cubic_spline(xlow, ylow),
        _cubic_spline(xup, yup),
    )

    xc = cosine_discretization(N)

    max_xt = min(
        maximum(xlow), maximum(xup)
    )
    isless_xt = @. xc <= max_xt
    xt = @view xc[isless_xt]

    Airfoil(
        _cubic_spline(
            xt, upper.(xt) .- lower.(xt)
        ),
        _cubic_spline(
            xc, (upper.(xc) .+ lower.(xc)) ./ 2
        ),
        LE_flap_chord, LE_flap_deflection,
        TE_flap_chord, TE_flap_deflection,
        forced_separation
    )

end

"""
```
    function angles(afl::Airfoil)
```

Get upper and lower angles of an airfoil surface with its camberline
"""
function angles(afl::Airfoil, η::AbstractVector)

    leflap_c = afl.LE_flap_chord
    teflap_c = afl.TE_flap_chord

    leflap_def = afl.LE_flap_deflection
    teflap_def = afl.TE_flap_deflection

    LE_flap_arm = @. clamp(
        η - leflap_c, - Inf, 0.0
    )
    TE_flap_arm = @. clamp(
        1.0 - teflap_c - η, - Inf, 0.0
    )

    csep = afl.forced_separation

    t = afl.thickness.(
        clamp.(η, 0.0, csep)
    )
    y = afl.camber.(η) .+ tand.(leflap_def) .* LE_flap_arm .+ tand.(teflap_def) .* TE_flap_arm

    u = @. y + t / 2
    l = @. y - t / 2

    @tullio θu[i + _] := atan((u[i + 1] - u[i]), (η[i + 1] - η[i]))
    @tullio θl[i + _] := atan((l[i + 1] - l[i]), (η[i + 1] - η[i]))

    (θu, θl)

end

using DelimitedFiles

"""
```
    Airfoil(fname::String; kwargs...) = Airfoil(
        readdlm(fname; skipstart = 1); kwargs...
    )
```

Read Selig format dat file
"""
Airfoil(fname::String; kwargs...) = Airfoil(
    readdlm(fname; skipstart = 1); kwargs...
)

"""
Dimensionalize x coordinates
"""
dimensionalize(xmin::Real, xmax::Real, x::AbstractVector) = (
    @. x * (xmax - xmin) + xmin
)

"""
Display an airfoil
"""
Base.show(io::IO, afl::Airfoil) = Base.show(io, "<Dreadnought.jl airfoil>")
