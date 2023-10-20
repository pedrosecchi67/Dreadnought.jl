function source_infl(
	cx,
	cy,
	cz,
	p1x,
	p1y,
	p1z,
	p2x,
	p2y,
	p2z;
	M∞::Real = 0.0
)

	β = sqrt(1.0 - M∞ ^ 2)

	ax = cx .- p1x'
	ay = (cy .- p1y') .* β
	az = (cz .- p1z') .* β
	
	bx = cx .- p2x'
	by = (cy .- p2y') .* β
	bz = (cz .- p2z') .* β

	na = @. sqrt(ax ^ 2 + ay ^ 2 + az ^ 2) + 1e-10
	nb = @. sqrt(bx ^ 2 + by ^ 2 + bz ^ 2) + 1e-10

    lx = @. ax - bx
    ly = @. ay - by
    lz = @. az - bz

    L = @. sqrt(lx ^ 2 + ly ^ 2 + lz ^ 2) + 1e-10

    lx = @. lx / L
    ly = @. ly / L
    lz = @. lz / L

    x1 = @. - (lx * ax + ly * ay + lz * az)
    x2 = @. - (lx * bx + ly * by + lz * bz)

    rx = @. ax + x1 * lx
    ry = @. ay + x1 * ly
    rz = @. az + x1 * lz

    R = @. sqrt(rx ^ 2 + ry ^ 2 + rz ^ 2) + 1e-10

    rx = @. rx / R
    ry = @. ry / R
    rz = @. rz / R

    vr = @. (
        x2 / (R * nb) - x1 / (R * na)
    ) / (4 * π)

    vl = @. (
        1.0 / nb - 1.0 / na
    ) / (4 * π)

    [
        (@. (vr * rx + vl * lx) / β ^ 2),
        (@. (vr * ry + vl * ly) / β),
        (@. (vr * rz + vl * lz) / β)
    ]

end
