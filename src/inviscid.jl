export inviscid_solution

include("forces.jl")

"""
```
    function inviscid_solution(
        acft::Aircraft;
        α::Real = 0.0,
        β::Real = 0.0,
        Mach::Real = 0.0,
        p::Real = 0.0,
        q::Real = 0.0,
        r::Real = 0.0,
        return_forces::Bool = true,
        momentum_reference = (0.0, 0.0, 0.0),
    )
```

Get inviscid solution for an aircraft and the forces acting upon it (two dictionaries in a tuple)
"""
function inviscid_solution(
    acft::Aircraft;
    α::Real = 0.0,
    β::Real = 0.0,
    Mach::Real = 0.0,
    p::Real = 0.0,
    q::Real = 0.0,
    r::Real = 0.0,
    return_forces::Bool = true,
    momentum_reference = (0.0, 0.0, 0.0),
)

    PGβ = sqrt(1.0 - Mach ^ 2)

    ω = [
        - p * 2 * acft.bref,
        q * 2 * acft.cref,
        - r * 2 * acft.bref,
    ]

    cx = acft.cx
    cy = acft.cy
    cz = acft.cz

    û = [
        cosd(α) * cosd(β),
        sind(β) * cosd(α),
        sind(α)
    ]

    u∞ = (
        û .+ let (dvx, dvy, dvz) = _cross(
            cx, cy, cz, ω...
        )
            permutedims(
                [
                    dvx dvy dvz
                ]
            )
        end
    )

    θu = acft.θ_upper
    θl = acft.θ_lower

    nxu = @. sin(θu)
    nxl = @. sin(θl)
    ncx = @. (nxu + nxl) / 2

    streamx = acft.ux
    streamy = acft.uy
    streamz = acft.uz
    
    nx = acft.nx
    ny = acft.ny
    nz = acft.nz

    uinfx, uinfy, uinfz = eachrow(u∞)

    uinf_stream = _dot(
        streamx, streamy, streamz,
        uinfx, uinfy, uinfz
    )

    λ = @. (nxu - nxl) * uinf_stream

    Δx = acft.S ./ acft.L

    Δsx = @. Δx / streamx

    σ = @. λ * Δsx

    RHS = _dot(
        nx, ny, nz, uinfx, uinfy, uinfz
    ) .- ncx .* uinf_stream

    Γ = - (acft.luAIC \ RHS) ./ PGβ

    Γ = Γ .* uinf_stream
    γ = @. Γ / Δx

    umsx = ((acft.source_AICx * σ) .+ uinfx) ./ PGβ
    umsy = ((acft.source_AICy * σ) .+ uinfy) ./ PGβ
    umsz = ((acft.source_AICz * σ) .+ uinfz) ./ PGβ

    ufx = ((acft.AICx * Γ) .+ uinfx)
    ufy = ((acft.AICy * Γ) .+ uinfy)
    ufz = ((acft.AICz * Γ) .+ uinfz)
    
    u_force = [
        ufx'; ufy'; ufz'
    ]
    ums = [umsx'; umsy'; umsz']

    fxi, fyi, fzi = vorticity_forces(acft, Γ, u_force)
    # fni = fxi .* nx .+ fyi .* ny .+ fzi .* nz
    # fci = fxi

    u_upper = @. umsx + γ * cos(θu) / 2
    u_lower = @. umsx - γ * cos(θl) / 2

    Cp_upper = @. 1.0 - u_upper ^ 2
    Cp_lower = @. 1.0 - u_lower ^ 2

    results = Dict{String, AbstractVecOrMat}()

    results["sigma"] = σ
    results["lambda"] = λ
    results["Gamma"] = Γ
    results["gamma"] = γ

    results["mean_velocity"] = ums
    results["force_velocity"] = u_force

    results["u_upper"] = u_upper
    results["u_lower"] = u_lower
    results["Cp_upper"] = Cp_upper
    results["Cp_lower"] = Cp_lower

    results["inviscid_force"] = [
        fxi'; fyi'; fzi'
    ]
    # results["inviscid_normal_force"] = fni
    # results["invicsid_chordwise_force"] = fci

    if return_forces
        forces = total_forces(
            acft, results;
            momentum_reference = momentum_reference,
            α = α, β = β, 
            p = p, q = q, r = r
        )

        return (results, forces)
    end

    results

end
