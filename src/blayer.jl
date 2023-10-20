using Tullio

module Eppler

    """
    Eppler's second model for semi-empirical transition to turbulence, as per 
    Eppler, R.: Turbulent Airfoils for General Aviation, Journal of Aircraft, Vol. 15, No. 2, 1978.
    """
    function lnReth(
        H, r::Real = 1.0
    )

        _h = @. clamp(H, 2.0, 4.0)
        _hst = @. 1.515 + 0.076 * (4.0 - _h) ^ 2 / _h

        @. 18.4 * _hst - 21.74 - 0.36 * r + 0.125 * (_hst - 1.573) ^ 2 - 0.5 / (4.1 - _h)

    end

end

module Thwaites

    function Cf(
        Λ, Reθ
    )

        _reth = @. clamp(Reθ, 10.0, 1e6)

        λ = @. clamp(Λ * _reth, -0.1, 0.15)

        Τ = @. 0.22 + 1.52 * λ - 5.0 * λ ^ 3 - 0.072 * λ ^ 2 / (λ + 0.18) ^ 2

        @. 2 * Τ / (Reθ + 1.0)

    end

    function H(
        Λ, Reθ
    )

        _reth = @. clamp(Reθ, 10.0, 1e6)

        λ = @. clamp(Λ * _reth, -0.1, 0.15)

        @. 2.61 - 4.1 * λ + 14.0 * λ ^ 3 + 0.56 * λ ^ 2 / (λ + 0.18) ^ 2

    end

end

module White

    function H(
        Λ, Reθ
    )

        _rth = @. clamp(Reθ, 1e2, 1e6)
        _Λ = @. clamp(Λ, -4.52859e-3, 4.5e-3)
        _L = @. log10(_rth)

        H = @. - 4.072 * log(_Λ + 4.5286e-3) / (- 0.1331 * _L ^ 2 + 1.3061 * _L + 6.0) - 1.085
        @. clamp(H, 1.0, 2.38)

    end

    function Cf(Λ, Reθ)

        _rth = @. clamp(Reθ, 1e2, 1e6)
        _L = @. log10(_rth)

        _H = @. H(Λ, Reθ)

        Cf = @. 0.3 * exp(- 1.33 * _H) / (
            _L ^ (1.74 + 0.31 * _H)
        )

    end

end

relu(x::Real) = x * (x > 0.0)

"""
Ponderation function to find 1 minus the ratio at which, in a segment of distance L,
property p variating from p1 to p2 reaches the value of zero.
"""
ponderation(
    f1::Real, f2::Real
) = (relu(f1) + relu(f2)) / (abs(f1) + abs(f2) + eps(typeof(f1)))

"""
Function to initialize the state vector and the intermediary variable vector of a BL calculation
"""
initial_state(
    ;
    transition::Bool = false,
) = (
    Real[
        0.0, # P
        (transition ? 1.0 : 0.0), # transition flag
        0.0 # separation flag
    ],
    Real[
        0.0, 0.0,
        2.5, 0.0,
        -1.0, -1.0,
        0.0
    ]
)

"""
Function to get intermediary values from a BL calculation
"""
function get_intermediary(
    q::AbstractVector, # state vector
    Re::Real, # Reynolds number
    ue::Real,
    due!dx::Real;
    M∞::Real = 0.0,
    r::Real = 0.0,
    H_separation::Real = 2.3,
    γ::Real = 1.4,
)

    qe = abs(ue) + eps(eltype(ue))

    M = M∞ * qe

    ρ = (
        (
            (1.0 + (γ - 1.0) * M∞ ^ 2 / 2)
        ) / 
        (
            (1.0 + (γ - 1.0) * M ^ 2 / 2)
        )
    ) ^ (1.0 / (γ - 1.0))

    P, trans, sep = q

    θ = P / (qe ^ 2 * ρ)

    Reθ = Re * θ * ρ * qe
    Reθ = clamp(Reθ, 10.0, 1e6)

    Λ = θ * due!dx / qe

    Hl, Cfl = Thwaites.H(Λ, Reθ), Thwaites.Cf(Λ, Reθ)
    Ht, Cft = White.H(Λ, Reθ), White.Cf(Λ, Reθ)

    H = Hl * (1.0 - trans) + Ht * trans
    Cf = Cfl * (1.0 - trans) + Cft * trans

    transition_threshold = log(Reθ) - Eppler.lnReth(Hl, r)
    separation_threshold = Ht - H_separation

    τ = Cf * qe ^ 2 * ρ / 2

    dU!dx = due!dx / qe
    
    y = [
        τ, dU!dx,
        H, Cf,
        transition_threshold, separation_threshold,
        θ
    ]

end

"""
March boundary layer properties across a mesh cell
"""
function march_cell(
    qim1::AbstractVector,
    yim1::AbstractVector,
    ueim1::Real,
    ue::Real,
    ds::Real;
    M∞::Real = 0.0,
    Re::Real = 1e6,
    r::Real = 0.0,
    H_separation::Real = 2.3,
    γ::Real = 1.4,
    n_iter::Int = 1,
    ω::Real = 1.0,
    squire_young::Bool = false, # uses squire-young correction for P
)

    due!dx = (ue - ueim1) / ds

    iter = y -> let (τ, dU!dx, H, Cf, _, _, θ) = y
        Pim1, transim1, sepim1 = qim1

        _, _, _, _, trans_th_im1, sep_th_im1 = yim1

        P = (
            dU!dx > 0.0 ?
            (Pim1 + τ * ds) / (1.0 + H * dU!dx * ds) :
            (Pim1 + τ * ds) * (1.0 - H * dU!dx * ds)
        )

        ynew = get_intermediary(
            [P, transim1, sepim1], Re, ue, due!dx;
            M∞ = M∞,
            r = r,
            H_separation = H_separation,
            γ = γ,
        )

        _, _, _, _, trans_th, sep_th, _ = ynew

        # max so as not to allow for relaminarization
        trans = max(
            transim1,
            ponderation(trans_th_im1, trans_th)
        )
        # sep = ponderation(sep_th_im1, sep_th)
        sep = max(
            sepim1, 
            ponderation(sep_th_im1, sep_th)
        )

        # only allow separation flag after transition
        sep = min(trans, sep)

        q = [P, trans, sep]

        (
            q, ynew
        )
    end

    q = copy(qim1)
    y = copy(yim1)

    for nit = 1:n_iter
        qnew, ynew = iter(y)

        @. q += (qnew - q) * ω
        @. y += (ynew - y) * ω
    end

    if squire_young
        P = q[1]
        θ = y[end]
        H = y[3]

        qe = abs(ue)

        P *= qe ^ ((H + 1.0) / 2)
        θ *= qe ^ ((H + 1.0) / 2)

        q[1] = P
        y[end] = θ
    end

    Pim1, _, _ = qim1
    P, _, _ = q

    drag = (P - Pim1) / ds

    (q, y, drag)

end

"""
March boundary layer solution around a one-dimensional surface with a 
single stagnation point. Selig (counter-clockwise) format is expected
"""
function march_around(
    ue::AbstractVector,
    ds::AbstractVector;
    M∞::Real = 0.0,
    Re::Real = 1e6,
    r::Real = 0.0,
    H_separation::Real = 2.0,
    γ::Real = 1.4,
    n_iter::Int = 1,
    ω::Real = 1.0,
    forced_transition::Bool = false,
    squire_young::Bool = true, # impose squire young correction at first and last cells
)

    istag = findfirst(
        i -> ue[i] >= 0.0,
        1:length(ue)
    )
    if isnothing(istag)
        istag = length(ue)
    end

    # ponderation factor to position the stagnation point exactly

    η = ponderation(
        ue[max(istag - 1, 1)],
        ue[istag]
    )

    dsattachment = (
        ds[max(istag - 1, 1)] +
        ds[istag]
    ) / 2

    dsneg = dsattachment * (1.0 - η) + ds[max(istag - 1, 1)] / 2
    dspos = dsattachment * η + ds[istag] / 2

    # negative

    qlast, ylast = initial_state(; transition = forced_transition,)

    tups_neg = [
        begin
            q, y, drag = march_cell(
                qlast, ylast,
                abs(ue[i + 1]), abs(ue[i]),
                (i == istag - 1 ? dsneg : ds[i]);
                M∞ = M∞, Re = Re,
                r = r, H_separation = H_separation,
                γ = γ, 
                n_iter = n_iter,
                ω = ω,
                squire_young = squire_young && (i == 1)
            )

            qlast .= q
            ylast .= y

            (q, y, drag)
        end for i = (istag - 1):-1:1
    ]

    # positive

    qlast, ylast = initial_state(; transition = forced_transition,)

    tups_pos = [
        begin
            q, y, drag = march_cell(
                qlast, ylast,
                (
                    i == 1 ?
                    sqrt(eps(eltype(ue))) :
                    abs(ue[i - 1])
                ), abs(ue[i]),
                (i == istag ? dspos : ds[i]);
                M∞ = M∞, Re = Re,
                r = r, H_separation = H_separation,
                γ = γ, 
                n_iter = n_iter,
                ω = ω,
                squire_young = squire_young && (i == length(ue))
            )

            qlast .= q
            ylast .= y

            (q, y, drag)
        end for i = istag:length(ue)
    ]

    tups = [
        tups_neg[end:-1:1]; tups_pos
    ]

    #=
    [
        0.0, # P
        (transition ? 1.0 : 0.0), # transition flag
        0.0 # separation flag
    ],
    y = [
        τ, dU!dx,
        H, Cf,
        transition_threshold, separation_threshold,
        θ
    ]
    =#
    results = mapreduce(
        t -> vcat(t...), hcat, tups
    )

    fetch_i = i -> results[i, :]

    θ = fetch_i(10)
    H = fetch_i(6)

    drag = fetch_i(11)

    δst = θ .* H

    (
        ue = copy(ue),
        P = fetch_i(1),
        τ = fetch_i(4),
        θ = θ,
        δst = δst,
        H = H,
        Cf = fetch_i(7),
        Cp = (@. 1.0 - ue ^ 2),
        drag = drag,
        transition_flag = fetch_i(2),
        separation_flag = fetch_i(3),
        transition_threshold = fetch_i(8),
        separation_threshold = fetch_i(9)
    )

end
