using DelimitedFiles

begin
    local soln, forces

    cd("4415_s20deg_wing")

    AR = 6.2
    Λ = 20.0
    Re = 2.1e6

    @info "Running NACA-4415, AR$AR, Λ = $Λ deg, Re $Re, Mach 0 test case"

    n4415 = Airfoil(
        NACA4(4415)...
    )

    wng_left = Surface(
        [
            - tand(Λ) * AR / 2, - AR / 2, 0.0
        ], 1.0,
        zeros(3), 1.0
    )
    wng_right = Surface(
        zeros(3), 1.0,
        [
            - tand(Λ) * AR / 2, AR / 2, 0.0
        ], 1.0
    )

    set_airfoils!(wng_left, n4415, n4415)
    set_airfoils!(wng_right, n4415, n4415)

    acft = Aircraft(wng_left, wng_right; Sref = AR, bref = AR, cref = 1.0)


    αs = collect(0.0:1.0:30.0)

    CLs = []
    CDs = []

    for α = αs
        soln, forces = viscous_solution(acft; α = α)

        total = forces["total"]

        push!(CLs, total["CL"])
        push!(CDs, total["CD"])
    end

    writedlm(
        "polar.dat",
        [
            ["AoA" "CL" "CD"];
            [
                αs CLs CDs
            ]
        ]
    )

    cd("..")
end
