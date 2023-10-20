export viscous_solution

include("blayer_coenen_drela.jl")

"""
```
    function viscous_solution(
        acft::Aircraft;
        α::Real = 0.0,
        β::Real = 0.0,
        Mach::Real = 0.0,
        p::Real = 0.0,
        q::Real = 0.0,
        r::Real = 0.0,
        Re::Real = 1e6,
        H_separation::Real = 2.0,
        γ::Real = 1.4,
        n_iter::Int = 1,
        ω::Real = 1.0,
        squire_young::Bool = true, # impose squire young correction at first and last cells
        momentum_reference = (0.0, 0.0, 0.0),
    )
```

Get viscous solution for an aircraft
"""
function viscous_solution(
    acft::Aircraft;
    α::Real = 0.0,
    β::Real = 0.0,
    Mach::Real = 0.0,
    p::Real = 0.0,
    q::Real = 0.0,
    r::Real = 0.0,
    Re::Real = 1e6,
    H_separation::Real = 2.0,
    γ::Real = 1.4,
    n_iter::Int = 1,
    ω::Real = 1.0,
    squire_young::Bool = true, # impose squire young correction at first and last cells
    momentum_reference = (0.0, 0.0, 0.0),
)

    Recorr = Re / acft.cref

    isoln, iforces = inviscid_solution(
        acft;
        α = α, β = β,
        p = p, q = q, r = r,
        Mach = Mach,
    )

    ds_upper = acft.S ./ acft.L ./ cos.(acft.θ_upper)
    ds_lower = acft.S ./ acft.L ./ cos.(acft.θ_lower)

    u_upper = isoln["u_upper"]
    u_lower = isoln["u_lower"]

    ϵ = sqrt(eps(eltype(u_upper)))

    u_upper = @. clamp(u_upper, ϵ, Inf)
    u_lower = @. clamp(u_lower, ϵ, Inf)

    translation_table = Dict(
        :θ => "theta",
        :δst => "deltastar",
        :τ => "tau"
    )

    bl_results_upper = []
    bl_results_lower = []

    for surf in acft.surfaces
        ds_u = ds_upper[surf]
        ds_l = ds_lower[surf]

        uu = u_upper[surf]
        ul = u_lower[surf]

        for (su, sl, ueu, uel) in zip(
            eachcol(ds_u), eachcol(ds_l),
            eachcol(uu), eachcol(ul)
        )
            soln_u = march_around(
                ueu, su;
                M∞ = Mach, Re = Recorr,
                Ncr = surf.Ncr,
                H_separation = H_separation,
                γ = γ,
                n_iter = n_iter, ω = ω,
                forced_transition = surf.forced_transition,
                squire_young = squire_young
            )

            soln_l = march_around(
                uel, sl;
                M∞ = Mach, Re = Recorr,
                Ncr = surf.Ncr,
                H_separation = H_separation,
                γ = γ,
                n_iter = n_iter, ω = ω,
                forced_transition = surf.forced_transition,
                squire_young = squire_young
            )

            push!(bl_results_upper, soln_u)
            push!(bl_results_lower, soln_l)
        end
    end

    results = isoln

    props = propertynames(first(bl_results_lower))

    for p in props
        if !haskey(translation_table, p)
            translation_table[p] = String(p)
        end
    end

    for p in props
        vu = mapreduce(
            blsoln -> getproperty(
                blsoln, p
            ),
            vcat,
            bl_results_upper
        )
        vl = mapreduce(
            blsoln -> getproperty(
                blsoln, p
            ),
            vcat,
            bl_results_lower
        )

        results[translation_table[p] * "_upper"] = vu 
        results[translation_table[p] * "_lower"] = vl
    end

    uf = results["force_velocity"]

    uf_norm = map(
        norm, eachcol(uf)
    )

    σ_profile = (results["drag_upper"] .* ds_upper .+ results["drag_lower"] .* ds_lower) ./ uf_norm

    fpx, fpy, fpz = source_forces(
        acft, σ_profile, uf
    )

    results["profile_drag_force"] = [
        fpx'; fpy'; fpz'
    ]

    forces = total_forces(
        acft, results;
        momentum_reference = momentum_reference,
        α = α, β = β, 
        p = p, q = q, r = r
    )

    for (k, v) in iforces["total"]
        forces["total"][k * "_inviscid"] = v
    end

    for (d, di) in zip(
        forces["sectional"],
        iforces["sectional"],
    )
        for (k, v) in di
            d[k * "_inviscid"] = v
        end
    end

    for (d, di) in zip(
        forces["surface"],
        iforces["surface"],
    )
        for (k, v) in di
            d[k * "_inviscid"] = v
        end
    end

    (results, forces)

end
