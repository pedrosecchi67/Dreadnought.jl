_disc_handle(v::AbstractVector) = v
_disc_handle(v::Int) = cosine_discretization(v)

export Surface

"""
```
    struct Surface
        chord_left::Real
        LE_left::AbstractVector

        chord_right::Real
        LE_right::AbstractVector

        ηc_left::AbstractVector # from 0 to 1, control points across the chord in leftmost section
        ηc_right::AbstractVector # from 0 to 1, control points across the chord in rightmost section
        ηb::AbstractVector # from 0 to 1, spacing of control points across the span

        p1x::AbstractMatrix # panel kin coordinates
        p1y::AbstractMatrix
        p1z::AbstractMatrix
        p2x::AbstractMatrix # panel kin coordinates
        p2y::AbstractMatrix
        p2z::AbstractMatrix
        p3x::AbstractMatrix # panel kin coordinates
        p3y::AbstractMatrix
        p3z::AbstractMatrix
        p4x::AbstractMatrix # panel kin coordinates
        p4y::AbstractMatrix
        p4z::AbstractMatrix

        pleft_x::AbstractMatrix
        pright_x::AbstractMatrix
        pleft_y::AbstractMatrix
        pright_y::AbstractMatrix
        pleft_z::AbstractMatrix
        pright_z::AbstractMatrix

        cx::AbstractMatrix
        cy::AbstractMatrix
        cz::AbstractMatrix

        support_x::AbstractMatrix
        support_y::AbstractMatrix
        support_z::AbstractMatrix

        ux::AbstractMatrix
        uy::AbstractMatrix
        uz::AbstractMatrix

        vx::AbstractMatrix
        vy::AbstractMatrix
        vz::AbstractMatrix

        nx::AbstractMatrix
        ny::AbstractMatrix
        nz::AbstractMatrix

        L::AbstractMatrix
        S::AbstractMatrix

        θ_upper::AbstractMatrix
        θ_lower::AbstractMatrix

        indices::AbstractMatrix

        forced_transition::Bool
        Ncr::Real

        quarter_chord_x::AbstractVector
        quarter_chord_y::AbstractVector
        quarter_chord_z::AbstractVector

        strip_width::AbstractVector
        chord::AbstractVector

        function Surface(
            LE_left::AbstractVector,
            chord_left::Real,
            LE_right::AbstractVector,
            chord_right::Real;
            ηc_left::Union{Int, AbstractVector} = 40,
            ηc_right::Union{Int, AbstractVector} = 40,
            ηb::Union{Int, AbstractVector} = 30, # if integers, translated to cosine discretizations
            forced_transition::Bool = false,
            Ncr::Real = 6.0,
        )

            #

        end
    end
```

Struct to define a wing surface
"""
struct Surface
    chord_left::Real
    LE_left::AbstractVector

    chord_right::Real
    LE_right::AbstractVector

    ηc_left::AbstractVector # from 0 to 1, control points across the chord in leftmost section
    ηc_right::AbstractVector # from 0 to 1, control points across the chord in rightmost section
    ηb::AbstractVector # from 0 to 1, spacing of control points across the span

    p1x::AbstractMatrix # panel kin coordinates
    p1y::AbstractMatrix
    p1z::AbstractMatrix
    p2x::AbstractMatrix # panel kin coordinates
    p2y::AbstractMatrix
    p2z::AbstractMatrix
    p3x::AbstractMatrix # panel kin coordinates
    p3y::AbstractMatrix
    p3z::AbstractMatrix
    p4x::AbstractMatrix # panel kin coordinates
    p4y::AbstractMatrix
    p4z::AbstractMatrix

    pleft_x::AbstractMatrix
    pright_x::AbstractMatrix
    pleft_y::AbstractMatrix
    pright_y::AbstractMatrix
    pleft_z::AbstractMatrix
    pright_z::AbstractMatrix

    cx::AbstractMatrix
    cy::AbstractMatrix
    cz::AbstractMatrix

    support_x::AbstractMatrix
    support_y::AbstractMatrix
    support_z::AbstractMatrix

    ux::AbstractMatrix
    uy::AbstractMatrix
    uz::AbstractMatrix

    vx::AbstractMatrix
    vy::AbstractMatrix
    vz::AbstractMatrix

    nx::AbstractMatrix
    ny::AbstractMatrix
    nz::AbstractMatrix

    L::AbstractMatrix
    S::AbstractMatrix

    θ_upper::AbstractMatrix
    θ_lower::AbstractMatrix

    indices::AbstractMatrix

    forced_transition::Bool
    Ncr::Real

    quarter_chord_x::AbstractVector
    quarter_chord_y::AbstractVector
    quarter_chord_z::AbstractVector

    strip_width::AbstractVector
    chord::AbstractVector

    function Surface(
        LE_left::AbstractVector,
        chord_left::Real,
        LE_right::AbstractVector,
        chord_right::Real;
        ηc_left::Union{Int, AbstractVector} = 40,
        ηc_right::Union{Int, AbstractVector} = 40,
        ηb::Union{Int, AbstractVector} = 30, # if integers, translated to cosine discretizations
        forced_transition::Bool = false,
        Ncr::Real = 6.0,
    )
    
        ηc_left = _disc_handle(ηc_left)
        ηc_right = _disc_handle(ηc_right)
        ηb = _disc_handle(ηb)

        @assert length(ηc_left) == length(ηc_right) "Incompatible discretization for tips at wing instance"

        LE = wing_interpolate(ηb, LE_left, LE_right)
        c = wing_interpolate(ηb, chord_left, chord_right)
        ηc = wing_interpolate(ηb, ηc_left, ηc_right)

        LEx = @view LE[1, :]
        LEy = @view LE[2, :]
        LEz = @view LE[3, :]

        px = LEx' .+ ηc .* c'
        py = repeat(
            LEy'; inner = (length(ηc_left), 1)
        )
        pz = repeat(
            LEz'; inner = (length(ηc_left), 1)
        )

        p1x = @view px[1:(end - 1), 1:(end - 1)]
        p1y = @view py[1:(end - 1), 1:(end - 1)]
        p1z = @view pz[1:(end - 1), 1:(end - 1)]
        p2x = @view px[1:(end - 1), 2:end]
        p2y = @view py[1:(end - 1), 2:end]
        p2z = @view pz[1:(end - 1), 2:end]
        p3x = @view px[2:end, 2:end]
        p3y = @view py[2:end, 2:end]
        p3z = @view pz[2:end, 2:end]
        p4x = @view px[2:end, 1:(end - 1)]
        p4y = @view py[2:end, 1:(end - 1)]
        p4z = @view pz[2:end, 1:(end - 1)]

        pleft_x = @. 0.25 * p4x + 0.75 * p1x
        pright_x = @. 0.25 * p3x + 0.75 * p2x
        pleft_y = @. 0.25 * p4y + 0.75 * p1y
        pright_y = @. 0.25 * p3y + 0.75 * p2y
        pleft_z = @. 0.25 * p4z + 0.75 * p1z
        pright_z = @. 0.25 * p3z + 0.75 * p2z

        paft_x = @. (p3x + p4x) / 2
        pfwd_x = @. (p1x + p2x) / 2
        paft_y = @. (p3y + p4y) / 2
        pfwd_y = @. (p1y + p2y) / 2
        paft_z = @. (p3z + p4z) / 2
        pfwd_z = @. (p1z + p2z) / 2

        cx = @. (pfwd_x * 0.25 + paft_x * 0.75)
        cy = @. (pfwd_y * 0.25 + paft_y * 0.75)
        cz = @. (pfwd_z * 0.25 + paft_z * 0.75)

        support_x = @. (pfwd_x * 0.75 + paft_x * 0.25)
        support_y = @. (pfwd_y * 0.75 + paft_y * 0.25)
        support_z = @. (pfwd_z * 0.75 + paft_z * 0.25)

        crx, cry, crz = begin
            crx = @. pright_x - pleft_x
            cry = @. pright_y - pleft_y
            crz = @. pright_z - pleft_z

            ncr = @. sqrt(crx ^ 2 + cry ^ 2 + crz ^ 2)

            (
                crx ./ ncr,
                cry ./ ncr,
                crz ./ ncr,
            )
        end

        nx, ny, nz = begin
            nx = zeros(size(crx))
            ny = - crz
            nz = copy(cry)

            nn = @. sqrt(nx ^ 2 + ny ^ 2 + nz ^ 2)

            (
                nx ./ nn,
                ny ./ nn,
                nz ./ nn,
            )
        end

        ux, uy, uz = _cross(
            crx, cry, crz,
            nx, ny, nz
        )

        # alias:
        vx = crx
        vy = cry
        vz = crz

        L = @. sqrt((pright_y - pleft_y) ^ 2 + (pright_z - pleft_z) ^ 2)
        S = @. L * (paft_x - pfwd_x)

        chord = sum.(
            eachcol(paft_x .- pfwd_x)
        )
        width = @view L[1, :]

        quarter_chord_x = chord ./ 4 .+ @view pfwd_x[1, :]
        quarter_chord_y = @view pfwd_y[1, :]
        quarter_chord_z = @view pfwd_z[1, :]

        new(
            chord_left, LE_left,
            chord_right, LE_right,
            ηc_left, ηc_right,
            ηb,
            p1x, p1y, p1z,
            p2x, p2y, p2z,
            p3x, p3y, p3z,
            p4x, p4y, p4z,
            pleft_x, pright_x,
            pleft_y, pright_y,
            pleft_z, pright_z,
            cx, cy, cz,
            support_x, support_y, support_z,
            ux, uy, uz,
            vx, vy, vz,
            nx, ny, nz,
            L, S,
            zeros(Real, size(cx)), zeros(Real, size(cx)),
            zeros(Int64, size(cx)),
            forced_transition, Ncr,
            quarter_chord_x, quarter_chord_y, quarter_chord_z,
            width, chord
        )

    end
end

export set_airfoils!

"""
Obtain average of values at adjacent panel kinks
"""
function _avg_panel_span(v::AbstractMatrix)

    (
        (@view v[:, 1:(end - 1)]) .+
        (@view v[:, 2:end])
    ) ./ 2

end

"""
```
    function set_airfoils!(
        surf::Surface, afl_left::Airfoil, afl_right::Airfoil;
        θ_left::Real = 0.0, θ_right::Real = 0.0, # incidences
    )
```

Set thickness and camber information according to airfoil geometries at the kinks of the wing
"""
function set_airfoils!(
    surf::Surface, afl_left::Airfoil, afl_right::Airfoil;
    θ_left::Real = 0.0, θ_right::Real = 0.0, # incidences
)

    incidence = deg2rad.(wing_interpolate(surf.ηb, θ_left, θ_right))

    θu_l, θl_l = angles(afl_left, surf.ηc_left)
    θu_r, θl_r = angles(afl_right, surf.ηc_right)

    θu = _avg_panel_span(wing_interpolate(surf.ηb, θu_l, θu_r) .+ incidence')
    θl = _avg_panel_span(wing_interpolate(surf.ηb, θl_l, θl_r) .+ incidence')

    surf.θ_upper .= θu
    surf.θ_lower .= θl

    surf

end

"""
Display a surface
"""
Base.show(io::IO, surf::Surface) = Base.show(io, "<Dreadnought.jl surface>")

"""
Interpolate a variable from its values at the leftmost and rightmost sections (vector)
"""
wing_interpolate(
    η::AbstractVector,
    vl::AbstractVector,
    vr::AbstractVector
) = vl .* (1.0 .- η)' .+ vr .* η'

"""
Interpolate a variable from its values at the leftmost and rightmost sections
"""
wing_interpolate(
    η::AbstractVector,
    vl::Real,
    vr::Real
) = vl .* (1.0 .- η) .+ vr .* η

"""
Interpolate a variable from its values at the leftmost and rightmost sections
"""
wing_interpolate(
    surf::Surface,
    vl,
    vr
) = wing_interpolate(surf.ηb, vl, vr)
