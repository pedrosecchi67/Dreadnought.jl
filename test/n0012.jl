begin
    @info "Running NACA-0012, AR5, Re 1e6, Mach 0 test case"

    n0012 = Airfoil(
        NACA4(0012)...
    )

    AR = 5.0

    wng = Surface(
        [
            0.0, - AR / 2, 0.0
        ], 1.0,
        [
            0.0, AR / 2, 0.0
        ], 1.0
    )

    set_airfoils!(wng, n0012, n0012)

    acft = Aircraft(wng; Sref = AR, bref = AR, cref = 1.0)


    soln, forces = viscous_solution(acft; α = 5.0)

    @show forces["total"]


    grid = vtk_grid("n0012", acft)

    vtk_record!(grid, soln)

    vtk_save(grid)
end
