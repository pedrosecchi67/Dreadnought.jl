rotate(
    pts::AbstractMatrix,
    θ::Real,
) = let M = [
    cosd(θ) (- sind(θ));
    sind(θ) cosd(θ)
]
    pts * M
end

"""
```
    function fowler_flap(
        airfoil::AbstractMatrix, # airfoil geometry
        flap::AbstractMatrix; # flap geometry (x from 0 to 1)
        flap_chord::Real = 0.3,
        gap::Real = 0.03,
        overlap::Real = 0.05,
        δ::Real = 35.0,
        remove::Bool = true, # remove lower surface of the airfoil where the slat is
    )
```

Generates fowler flap geometry.
Returns new airfoil and fowler flap geometries in Selig format, respectively 
"""
function fowler_flap(
    airfoil::AbstractMatrix,
    flap::AbstractMatrix;
    flap_chord::Real = 0.3,
    gap::Real = 0.03,
    overlap::Real = 0.05,
    δ::Real = 35.0,
    remove::Bool = true,
)

    im = size(airfoil, 1)
    if remove
        for i = size(airfoil, 1):-1:1
            if airfoil[i, 1] < airfoil[1, 1] - flap_chord
                im = i
                break
            end
        end
    end

    airfoil = airfoil[1:im, :]

    x0 = (
        remove ?
        airfoil[1, 1] : airfoil[end, 1]
    )
    y0 = (
        remove ?
        airfoil[1, 2] : airfoil[end, 2]
    )

    flap = rotate(flap, δ) .* flap_chord

    xmin = minimum(flap[:, 1])
    ymax = maximum(flap[:, 2])

    flap[:, 1] .+= x0 - overlap - xmin
    flap[:, 2] .+= y0 - gap - ymax

    (airfoil, flap)

end

"""
```
    function slot(
        airfoil::AbstractMatrix; # original airfoil in Selig format
        c_upper::Real = 0.2, # x position of slot start, on the upper surface
        c_lower::Real = 0.1, # x position of slot start, on the lower surface
        gap::Real = 0.05,
        η::Real = 0.5, # ponderation factor to estabilish the WUSS (Wing Under Slat Surface)
    ) # between a line connecting `c_lower`, `c_upper` points and the leading edge curve
```

Generate slot geometry.
Returns new slot and airfoil points, respectively, in Selig format
"""
function slot(
    airfoil::AbstractMatrix;
    c_upper::Real = 0.2,
    c_lower::Real = 0.1,
    gap::Real = 0.05,
    η::Real = 0.5,
)

    iu = 1
    for i = 1:size(airfoil, 1)
        if airfoil[i, 1] < c_upper
            iu = i

            break
        end
    end

    il = size(airfoil, 1)
    for i = size(airfoil, 1):(-1):1
        if airfoil[i, 1] < c_lower
            il = i

            break
        end
    end

    xju = airfoil[iu, 1]
    yju = airfoil[iu, 2]

    xjl = airfoil[il, 1]
    yjl = airfoil[il, 2]

    slot_geom = airfoil[iu:il, :]

    m = [0.0]
    for i = 2:size(slot_geom, 1)
        push!(
            m,
            m[end] + norm(slot_geom[i, :] .- slot_geom[i - 1, :])
        )
    end

    m ./= m[end]

    surf = [
        (xju .* (1.0 .- m) .+ xjl .* m) (yju .* (1.0 .- m) .+ yjl .* m)
    ]

    airfoil[iu:il, :] .= surf .* (1.0 - η) .+ slot_geom .* η

    slot_geom[:, 1] .-= gap

    (slot_geom, airfoil)

end

"""
```
    function slat(
        airfoil::AbstractMatrix; # original airfoil in Selig format
        c_upper::Real = 0.2, # x position of slot start, on the upper surface
        c_lower::Real = 0.05, # x position of slot start, on the lower surface
        gap::Real = 0.025,
        overlap::Real = 0.015,
        δ::Real = 35.0,
        η::Real = 0.5, # ponderation factor to estabilish the WUSS (Wing Under Slat Surface)
    ) # between a line connecting `c_lower`, `c_upper` points and the leading edge curve
```

Generate slat geometry.
Returns new slat and airfoil points, respectively, in Selig format
"""
function slat(
    airfoil::AbstractMatrix;
    c_upper::Real = 0.2,
    c_lower::Real = 0.05,
    gap::Real = 0.025,
    overlap::Real = 0.015,
    δ::Real = 35.0,
    η::Real = 0.5,
)

    airfoil = copy(airfoil)

    iu = 1
    for i = 1:size(airfoil, 1)
        if airfoil[i, 1] < c_upper
            iu = i

            break
        end
    end

    il = size(airfoil, 1)
    for i = size(airfoil, 1):(-1):1
        if airfoil[i, 1] < c_lower
            il = i

            break
        end
    end

    xju = airfoil[iu, 1]
    yju = airfoil[iu, 2]

    xjl = airfoil[il, 1]
    yjl = airfoil[il, 2]

    slot_geom = airfoil[iu:il, :]

    m = [0.0]
    for i = 2:size(slot_geom, 1)
        push!(
            m,
            m[end] + norm(slot_geom[i, :] .- slot_geom[i - 1, :])
        )
    end

    m ./= m[end]

    surf = [
        (xju .* (1.0 .- m) .+ xjl .* m) (yju .* (1.0 .- m) .+ yjl .* m)
    ]

    airfoil[iu:il, :] .= surf .* (1.0 - η) .+ slot_geom .* η

    slot_geom = rotate(slot_geom, - δ)

    ile = argmin(airfoil[:, 1])
    xle = airfoil[ile, 1]
    yle = airfoil[ile, 2]

    xn = xle - overlap
    yn = yle + sqrt(gap ^ 2 - overlap ^ 2)

    slot_geom[:, 1] .+= (xn - slot_geom[1, 1])
    slot_geom[:, 2] .+= (yn - slot_geom[1, 2])

    (slot_geom, airfoil)

end

"""
```
    function spoiler(
        pts::AbstractMatrix;
        δ::Real = 70.0,
        spoiler_origin::Real = 0.7,
        spoiler_chord::Real = 0.25,
    )
```

Add spoiler to an airfoil geometry and return the result
"""
function spoiler(
    pts::AbstractMatrix;
    δ::Real = 70.0,
    spoiler_origin::Real = 0.7,
    spoiler_chord::Real = 0.25,
)

    pts = copy(pts)

    for i = 1:size(pts, 1)
        if pts[i, 1] < spoiler_origin + spoiler_chord
            pts = pts[i:end, :]
            
            break
        end
    end

    root = spoiler_origin

    iu = 0
    for i = 1:size(pts, 1)
        if pts[i, 1] < root
            iu = i
            
            break
        end
    end

    #=
    il = 0
    for i = size(pts, 1):(-1):1
        if pts[i, 1] < root
            il = i
            
            break
        end
    end
    =#

    pts[1:iu, :] .= let p0 = pts[iu, :]
        p0' .+ rotate(pts[1:iu, :] .- p0', - δ)
    end

    pts

end

"""
```
    function split_flap(
        pts::AbstractMatrix;
        δ::Real = 40.0,
        flap_origin::Real = 0.75,
        flap_chord::Real = 0.25,
    )
```

Add split flap to the airfoil geometry and return the altered geom.
"""
function split_flap(
    pts::AbstractMatrix;
    δ::Real = 40.0,
    flap_origin::Real = 0.75,
    flap_chord::Real = 0.25,
)

    pts = copy(pts)

    for i = size(pts, 1):(-1):1
        if pts[i, 1] < flap_origin + flap_chord
            pts = pts[1:i, :]
            
            break
        end
    end

    root = flap_origin

    #=
    iu = 0
    for i = 1:size(pts, 1)
        if pts[i, 1] < root
            iu = i
            
            break
        end
    end
    =#

    il = 0
    for i = size(pts, 1):(-1):1
        if pts[i, 1] < root
            il = i
            
            break
        end
    end

    pts[il:end, :] .= let p0 = pts[il, :]
        p0' .+ rotate(pts[il:end, :] .- p0', δ)
    end

    pts

end


"""
```
    function plain_flap(
        pts::AbstractMatrix;
        δ::Real = 30.0,
        flap_chord::Real = 0.3,
    )
```

Add plain flap to airfoil geometry and return the deformed geom.
"""
function plain_flap(
    pts::AbstractMatrix;
    δ::Real = 30.0,
    flap_chord::Real = 0.3,
)

    pts = copy(pts)

    root = pts[1, 1] - flap_chord

    iu = 0
    for i = 1:size(pts, 1)
        if pts[i + 1, 1] < root
            iu = i
            
            break
        end
    end

    il = 0
    for i = size(pts, 1):(-1):1
        if pts[i - 1, 1] < root
            il = i
            
            break
        end
    end

    pts[1:iu, 2] .-= sind(δ) .* (pts[1:iu, 1] .- root)
    pts[il:end, 2] .-= sind(δ) .* (pts[il:end, 1] .- root)

    pts[1:iu, 1] .-= (1.0 - cosd(δ)) .* (pts[1:iu, 1] .- root)
    pts[il:end, 1] .-= (1.0 - cosd(δ)) .* (pts[il:end, 1] .- root)

    pts

end


"""
```
    function leading_edge_flap(
        pts::AbstractMatrix;
        δ::Real = 20.0,
        flap_chord::Real = 0.2,
    )
```

Add leading edge flap to an airfoil and return the deformed geometry
"""
function leading_edge_flap(
    pts::AbstractMatrix;
    δ::Real = 20.0,
    flap_chord::Real = 0.2,
)

    pts = copy(pts)

    ile = argmin(pts[:, 1])
    root = pts[ile, 1] + flap_chord

    iu = 0
    for i = 1:size(pts, 1)
        if pts[i, 1] < root
            iu = i
            
            break
        end
    end

    il = 0
    for i = size(pts, 1):(-1):1
        if pts[i, 1] < root
            il = i
            
            break
        end
    end

    pts[iu:il, 2] .+= sind(δ) .* (pts[iu:il, 1] .- root)

    pts[iu:il, 1] .-= (1.0 - cosd(δ)) .* (pts[iu:il, 1] .- root)

    pts

end
