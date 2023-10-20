using DelimitedFiles

begin
    cd("4412_airfoil")

    AR = 100.0
    Re = 9e6

    @info "Running NACA-4412, AR$AR, Re $Re, Mach 0 test case"

    n4412 = Airfoil(
        NACA4(4412)...
    )

    wng = Surface(
        [
            -0.25, - AR / 2, 0.0
        ], 1.0,
        [
            -0.25, AR / 2, 0.0
        ], 1.0
    )

    set_airfoils!(wng, n4412, n4412)

    acft = Aircraft(wng; Sref = AR, bref = AR, cref = 1.0)


    soln, forces = viscous_solution(acft; α = 15.0)

    @show forces["total"]


    grid = vtk_grid("n4412", acft)

    vtk_record!(grid, soln)

    vtk_save(grid)

    αs = collect(0.0:1.0:25.0)

    CLs = []
    CDs = []
    Cms = []

    for α = αs
        local soln, forces

        soln, forces = viscous_solution(acft; α = α)

        total = forces["total"]

        push!(CLs, total["CL"])
        push!(CDs, total["CD"] - total["CD_inviscid"])
        push!(Cms, total["Cm"])
    end

    writedlm(
        "polar.dat",
        [
            ["AoA" "CL" "CD" "Cm"];
            [
                αs CLs CDs Cms
            ]
        ]
    )

    cd("..")
end
