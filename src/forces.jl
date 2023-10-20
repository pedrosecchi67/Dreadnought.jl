"""
```
    function vorticity_forces(
        acft::Aircraft, Γ::AbstractVector, us::AbstractMatrix
    )
```

Obtain forces acting at nodes given a vorticity distribution and 
flow velocities at bound vortices
"""
function vorticity_forces(
    acft::Aircraft, Γ::AbstractVector, us::AbstractMatrix
)

    p1x = acft.pleft_x
    p1y = acft.pleft_y
    p1z = acft.pleft_z
    p2x = acft.pright_x
    p2y = acft.pright_y
    p2z = acft.pright_z

    fx, fy, fz = _cross(
        eachrow(us)...,
        Γ .* (p2x .- p1x),
        Γ .* (p2y .- p1y),
        Γ .* (p2z .- p1z),
    )

    (fx, fy, fz)

end

"""
```
    function source_forces(
        acft::Aircraft, σ::AbstractVector, us::AbstractMatrix
    )
```

Obtain forces acting at nodes given a source distribution and 
flow velocities at bound vortices
"""
function source_forces(
    acft::Aircraft, σ::AbstractVector, us::AbstractMatrix
)

    L = acft.L

    map(
        u -> (L .* u .* σ),
        eachrow(us)
    )

end

"""
```
    TE_thickness(surf::Surface)
```

Get vector with spanwise trailing edge thicknesses
"""
TE_thickness(surf::Surface) = vec(
    sum(
        (tan.(surf.θ_upper) .- tan.(surf.θ_lower)) .* surf.S ./ surf.L; dims = 1
    )
)

#=
"""
Obtain lift rescaling factor from local separation flags and Cp distribution
"""
function kirchhoff_factor(
    dS::AbstractVector,
    separation_flag::AbstractVector,
    Cp_upper::AbstractVector,
    Cp_lower::AbstractVector,
)

    dCp = Cp_upper .- Cp_lower

    K = abs(
        sum(
            dCp .* dS .* (1.0 .- separation_flag)
        )
    ) / (
        abs(
            sum(dCp .* dS)
        ) + sqrt(eps(eltype(dCp)))
    )

    K * 0.75 + 0.25

end
=#

"""
Obtain variation in sectional lift and drag due to the action of a decambering flap upon separated regions
"""
function decambering_flap(
    surf::Surface,
    upper_flag::AbstractMatrix,
    lower_flag::AbstractMatrix,
    α_local::AbstractVector
)

    dS = surf.S

    θu = - surf.θ_upper .+ deg2rad.(α_local)'
    θl = - surf.θ_lower .+ deg2rad.(α_local)'
    
    dCL = - vec(
        sum(
            upper_flag .* θu .* dS; dims = 1,
        ) + sum(
            lower_flag .* θl .* dS; dims = 1,
        )
    ) .* π ./ vec(sum(dS; dims = 1))

    dCL

end

"""
```
    function surface_forces(
        acft::Aircraft,
        surf::Surface,
        soln::AbstractDict;
        α::Real = 0.0,
        β::Real = 0.0,
        p::Real = 0.0,
        q::Real = 0.0,
        r::Real = 0.0,
        momentum_reference = (0.0, 0.0, 0.0),
    )
```

Obtain sectional and total forces for a lifting surface (returned as tuple of dicts)
"""
function surface_forces(
    acft::Aircraft,
    surf::Surface,
    soln::AbstractDict;
    α::Real = 0.0,
    β::Real = 0.0,
    p::Real = 0.0,
    q::Real = 0.0,
    r::Real = 0.0,
    momentum_reference = (0.0, 0.0, 0.0),
)

    results = Dict{String, AbstractVector}()

    getprop = (p; fline = true, dim = 0) -> (
        let v = let gp = getproperty(acft, p)
            (
                dim == 0 ?
                gp[surf] : gp[dim, :][surf]
            )
        end
            (
                fline ?
                v[1, :] : v
            )
        end
    )
    getsoln = (v; dim = 0) -> let x = soln[v]
        (
            dim == 0 ?
            x[surf] : x[dim, :][surf]
        )
    end

    ny = getprop(:ny)
    nz = getprop(:nz)

    qcx = surf.quarter_chord_x
    qcy = surf.quarter_chord_y
    qcz = surf.quarter_chord_z

    results["x"] = qcx
    results["y"] = qcy
    results["z"] = qcz

    c = surf.chord
    w = surf.strip_width

    results["chord"] = c

    is_viscous = haskey(soln, "P_upper")

    fx = getsoln("inviscid_force"; dim = 1) .+ (
        is_viscous ? getsoln("profile_drag_force"; dim = 1) : 0.0
    )
    fy = getsoln("inviscid_force"; dim = 2) .+ (
        is_viscous ? getsoln("profile_drag_force"; dim = 2) : 0.0
    )
    fz = getsoln("inviscid_force"; dim = 3) .+ (
        is_viscous ? getsoln("profile_drag_force"; dim = 3) : 0.0
    )

    bx = getprop(:support_x; fline = false)

    my = (fy .* ny' .+ fz .* nz') .* (qcx' .- bx)

    Fx, Fy, Fz, My = map(
        v -> vec(sum(v; dims = 1)),
        (fx, fy, fz, my)
    )

    sx, sy, sz = _cross(
        0.0, ny, nz, 1.0, 0.0, 0.0
    )

    ω = [
        - p * 2 * acft.bref,
        q * 2 * acft.cref,
        - r * 2 * acft.bref,
    ]

    û = [
        cosd(α) * cosd(β),
        sind(β) * cosd(α),
        sind(α)
    ]

    uinfx, uinfy, uinfz = map(
        (rot, uh) -> rot .+ uh,
        _cross(
            qcx, qcy, qcz, ω...
        ),
        û
    )

    α_local = atand.(
        uinfy .* ny .+ uinfz .* nz,
        uinfx
    )
    
    Q = @. sqrt(uinfx ^ 2 + uinfy ^ 2 + uinfz ^ 2) / 2
    dS = w .* c

    N = Fy .* ny .+ Fz .* nz
    C = Fx
    S = sx .* Fx .+ sy .* Fy .+ sz .* Fz

    Cn = N ./ (Q .* dS)
    Cc = C ./ (Q .* dS)
    Cs = S ./ (Q .* dS)
    Cm = My ./ (Q .* dS .* c)

    Cl = Cn .* cosd.(α_local) .- Cc .* sind.(α_local)
    Cd = Cn .* sind.(α_local) .+ cosd.(α_local) .* Cc

    if is_viscous
        Cd_profile = let (Pu, Pl) = (getsoln("P_upper"), getsoln("P_lower"))
            PTE = Pu[end, :] .+ Pl[end, :]

            2 .* PTE ./ Q ./ c
        end

        results["Cd_profile"] = Cd_profile

        uTE = sqrt.(
            getsoln("mean_velocity"; dim = 1)[end, :] .^ 2 .+
            getsoln("mean_velocity"; dim = 2)[end, :] .^ 2 .+
            getsoln("mean_velocity"; dim = 3)[end, :] .^ 2
        )
        tTE = TE_thickness(surf)

        #=
        Cp_upper = getsoln("Cp_upper")
        Cp_lower = getsoln("Cp_lower")
        =#

        panel_S = getprop(:S; fline = false)

        f_from_flag = (pS, sf) -> 1.0 - sum(pS .* sf) / sum(pS)

        upper_flag = getsoln("separation_flag_upper")
        lower_flag = getsoln("separation_flag_lower")
        flag = max.(upper_flag, lower_flag)

        f = map(
            f_from_flag,
            eachcol(panel_S), eachcol(flag)
        )

        results["f"] = f

        #=
        K = map(
            kirchhoff_factor,
            eachcol(panel_S),
            eachcol(flag),
            eachcol(Cp_upper), eachcol(Cp_lower)
        )
        # K = @. (sqrt(f)  + 1.0) ^ 2 / 4

        Cl_new = @. Cl * K

        dCl = Cl_new .- Cl
        =#

        dCl = decambering_flap(
            surf, upper_flag, lower_flag, α_local
        )

        Cl = Cl .+ dCl

        Cm = Cm .- dCl .* (1.0 .- f ./ 2 .- 0.25)

        Cd = Cd .+ tTE .* uTE .^ 2 .- Cd_profile .* (1.0 .- f) .+ sind.(α_local) .^ 2 .* π .* (1.0 .- f) ./ 2

        # Cl = Cl_new
    end

    Cn = Cl .* cosd.(α_local) .+ Cd .* sind.(α_local)
    Cc = - Cl .* sind.(α_local) .+ cosd.(α_local) .* Cd

    results["Cl"] = Cl
    results["Cd"] = Cd

    results["Cm"] = Cm
    results["alpha_local"] = α_local

    N = Cn .* Q .* dS
    C = Cc .* Q .* dS
    S = Cs .* Q .* dS
    My = Cm .* Q .* dS .* c

    Fx = C .+ S .* sx
    Fy = N .* ny .+ S .* sy
    Fz = N .* nz .+ S .* sz

    mrx, mry, mrz = momentum_reference

    Mx, dMy, Mz = _cross(
        qcx .- mrx, qcy .- mry, qcz .- mrz,
        Fx, Fy, Fz
    )

    My = My .+ dMy

    Fx, Fy, Fz, Mx, My, Mz = sum.(
        (
            Fx, Fy, Fz, Mx, My, Mz
        )
    )

    Q = 0.5

    S = acft.Sref
    b = acft.bref
    c = acft.cref

    F = [Fx, Fy, Fz]

    D = F ⋅ û
    L = norm(F .- D .* û)

    CX, CY, CZ = (
        - Fx / Q / S,
        Fy / Q / S,
        - Fz / Q / S,
    )

    CD = D / Q / S
    CL = L / Q / S
    Cl = Mx / Q / S / b
    Cm = My / Q / S / c
    Cn = Mz / Q / S / b

    total = Dict(
        "CX" => CX,
        "CY" => CY,
        "CZ" => CZ,
        "CD" => CD,
        "CL" => CL,
        "Cl" => Cl,
        "Cm" => Cm,
        "Cn" => Cn,
    )

    (results, total)

end

"""
```
    function total_forces(
        acft::Aircraft, soln::AbstractDict;
        kwargs...
    )
```

Obtain total forces acting upon an aircraft
"""
function total_forces(
    acft::Aircraft, soln::AbstractDict;
    kwargs...
)

    total = Dict{String, Real}()
    sectional = AbstractDict[]
    surf_forces = AbstractDict[]

    for surf in acft.surfaces
        sf, tf = surface_forces(
            acft, surf, soln; kwargs...
        )

        push!(sectional, sf)
        push!(surf_forces, tf)

        for (k, v) in tf
            if !haskey(total, k)
                total[k] = 0.0
            end
            
            total[k] += v
        end
    end

    Dict(
        "sectional" => sectional,
        "total" => total,
        "surface" => surf_forces,
    )

end
