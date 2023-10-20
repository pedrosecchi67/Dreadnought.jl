function _cross(
    ux, uy, uz,
    vx, vy, vz
)

    [
        (@. uy * vz - uz * vy),
        (@. uz * vx - ux * vz),
        (@. ux * vy - uy * vx),
    ]

end

function _dot(
    ux, uy, uz,
    vx, vy, vz
)

    @. ux * vx + uy * vy + uz * vz

end

function _norm(ux, uy, uz)

    sqrt.(ux .^ 2 .+ uy .^ 2 .+ uz .^ 2)

end

"""
Flatten an array to one dimension
"""
flatten(a::AbstractArray) = reshape(a, (length(a),))
