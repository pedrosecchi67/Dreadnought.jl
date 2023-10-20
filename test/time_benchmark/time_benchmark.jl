using DelimitedFiles

direct = (Nc, Nb, npre = 2, ntests = 4) -> begin
    tm = 0.0

    for i = 1:(npre + ntests)
        t0 = time()

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
            ], 1.0;
            ηc_left = Nc,
            ηc_right = Nc,
            ηb = Nb
        )

        set_airfoils!(wng, n0012, n0012)

        acft = Aircraft(wng; Sref = AR, bref = AR, cref = 1.0)


        soln = viscous_solution(acft; α = 5.0)

        if i > npre
            tm += time() - t0
        end
    end

    tm / ntests
end

precomputed = (Nc, Nb, npre = 2, ntests = 4) -> begin
    tm = 0.0

    wng = Surface(
        [
            0.0, - AR / 2, 0.0
        ], 1.0,
        [
            0.0, AR / 2, 0.0
        ], 1.0;
        ηc_left = Nc,
        ηc_right = Nc,
        ηb = Nb
    )

    acft = Aircraft(wng; Sref = AR, bref = AR, cref = 1.0)

    for i = 1:(npre + ntests)
        t0 = time()

        n0012 = Airfoil(
            NACA4(0012)...
        )

        AR = 5.0

        set_airfoils!(wng, n0012, n0012)


        soln = viscous_solution(acft; α = 5.0)

        if i > npre
            tm += time() - t0
        end
    end

    tm / ntests
end

begin
    cd("time_benchmark/")

    cases = [
        (10, 10),
        (20, 20),
        (30, 30),
        (40, 40),
        (50, 50),
        (60, 60),
    ]

    results = Any[
        ["Nc", "Nb", "N", "Direct(s)", "Precomputed(s)"]
    ]

    for case in cases
        Nc, Nb = case

        @info "Testing case Nc = $Nc, Nb = $Nb"

        td = direct(Nc, Nb)
        tp = precomputed(Nc, Nb)

        @info "Direct = $td s, Precomputed = $tp s"

        push!(
            results, [Nc, Nb, Nc * Nb, td, tp]
        )
    end

    results = permutedims(hcat(results...))

    writedlm("results.dat", results)

    cd("..")
end