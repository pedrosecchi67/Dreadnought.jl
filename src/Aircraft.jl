export Aircraft

"""
```
    struct Aircraft
        surfaces::Vector{Surface}

        Sref::Real
        bref::Real
        cref::Real

        Aircraft(
            surfaces::Surface...;
            Sref::Real = 1.0,
            cref::Real = 1.0,
            bref::Real = 1.0,
        ) = new(
            collect(surfaces),
            Sref, bref, cref
        )
    end
```

Struct to model an aircraft
"""
mutable struct Aircraft
    surfaces::Vector{Surface}

    Sref::Real
    bref::Real
    cref::Real

    source_AICx::Union{AbstractMatrix, Nothing}
    source_AICy::Union{AbstractMatrix, Nothing}
    source_AICz::Union{AbstractMatrix, Nothing}

    AICx::Union{AbstractMatrix, Nothing}
    AICy::Union{AbstractMatrix, Nothing}
    AICz::Union{AbstractMatrix, Nothing}

    luAIC::Union{LU, Nothing}

    function Aircraft(
        surfaces::Surface...;
        Sref::Real = 1.0,
        cref::Real = 1.0,
        bref::Real = 1.0,
        AIC_calculate::Bool = true,
    )

        n0 = 0
        for surf in surfaces
            surf.indices .= reshape(
                collect((n0 + 1):(n0 + length(surf.indices))),
                size(surf.indices)
            )

            n0 += length(surf.indices)
        end
    
        n = new(
            collect(surfaces),
            Sref, bref, cref,
            nothing, nothing, nothing,
            nothing, nothing, nothing,
            nothing,
        )

        if AIC_calculate
            AIC_calculate!(n)
        end

        n

    end

end

include("source_kernel.jl")
include("vlm_kernel.jl")

"""
```
    function AIC_calculate!(acft::Aircraft)
```

Calculate aerodynamic influence coefficient matrix
"""
function AIC_calculate!(acft::Aircraft)

    p1x = acft.pleft_x
    p1y = acft.pleft_y
    p1z = acft.pleft_z

    p2x = acft.pright_x
    p2y = acft.pright_y
    p2z = acft.pright_z

    cx = acft.cx
    cy = acft.cy
    cz = acft.cz

    nx = acft.nx
    ny = acft.ny
    nz = acft.nz

    bx = acft.support_x
    by = acft.support_y
    bz = acft.support_z

    # sources
    
    source_AICx, source_AICy, source_AICz = source_infl(
        cx, cy, cz,
        p1x, p1y, p1z,
        p2x, p2y, p2z,
    )

    acft.source_AICx = source_AICx
    acft.source_AICy = source_AICy
    acft.source_AICz = source_AICz

    # inverse AIC

    Ax, Ay, Az = hshoe_kernel(
        cx, cy, cz, 1.0, 0.0, 0.0, p1x, p1y, p1z, p2x, p2y, p2z,
    )

    acft.luAIC = lu(
        Ax .* nx .+ Ay .* ny .+ Az .* nz
    )

    # AIC

    Ax, Ay, Az = hshoe_kernel(
        bx, by, bz, 1.0, 0.0, 0.0, p1x, p1y, p1z, p2x, p2y, p2z,
    )

    acft.AICx = Ax
    acft.AICy = Ay
    acft.AICz = Az

end

"""
```
    function AIC_wipe!(acft::Aircraft)
```

Wipe and deallocate aerodynamic influence coefficient matrix
"""
function AIC_wipe!(acft::Aircraft)

    acft.AICx = nothing
    acft.AICy = nothing
    acft.AICz = nothing

    acft.source_AICx = nothing
    acft.source_AICy = nothing
    acft.source_AICz = nothing

end

"""
Concatenate the matrix-shaped values of properties across aerodynamic surfaces
"""
surfcat(vals::AbstractMatrix...) = mapreduce(
    vec, vcat, vals
)

"""
Get index, and concatenate from surfaces if absent
"""
Base.getproperty(
    acft::Aircraft, s::Symbol
) = (
    hasfield(Aircraft, s) ?
    getfield(acft, s) :
    surfcat(
        map(
            surf -> getproperty(surf, s),
            acft.surfaces
        )...
    )
)

"""
Get values of a property at a given aerodynamic surface
"""
Base.getindex(
    v::AbstractVector, surf::Surface
) = v[surf.indices]

using WriteVTK

export vtk_grid, vtk_save

"""
Mix a set of vectors
"""
_mixvs(vs::AbstractVector...) = vec(
    mapreduce(
        x -> x', vcat, vs
    )
)

"""
```
    function vtk_grid(fname::String, acft::Aircraft)
```

Create VTK grid file using WriteVTK.jl
"""
function vtk_grid(fname::String, acft::Aircraft)

    nc = mapreduce(
        surf -> length(surf.indices),
        +,
        acft.surfaces
    )

    pts = mapreduce(
        permutedims, vcat,
        (
            _mixvs(acft.p1x, acft.p2x, acft.p3x, acft.p4x),
            _mixvs(acft.p1y, acft.p2y, acft.p3y, acft.p4y),
            _mixvs(acft.p1z, acft.p2z, acft.p3z, acft.p4z),
        )
    )
    cells = reshape(collect(1:(4 * nc)), (4, nc))

    cells = map(
        c -> MeshCell(
            VTKCellTypes.VTK_QUAD, c
        ),
        eachcol(cells)
    )

    WriteVTK.vtk_grid(
        fname, pts, cells
    )

end

export vtk_record!

"""
```
    function vtk_record!(
        grid, results::AbstractDicts;
        suffix::String = "", # added after dict. keys
    )
```

Record dictionary of results in VTK grid
"""
function vtk_record!(
    grid, results::AbstractDict;
    suffix::String = "", # added after dict. keys
)

    for (k, v) in results
        grid[k * suffix] = Float64.(v)
    end

    grid

end
