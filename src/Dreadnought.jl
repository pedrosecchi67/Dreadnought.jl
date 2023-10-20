module Dreadnought

    using Interpolations
    using LinearAlgebra
    using Tullio

    include("vector_ops.jl")

    include("Airfoil.jl")
    include("Surface.jl")
    include("Aircraft.jl")

    include("inviscid.jl")
    include("viscous.jl")

    # don't forget to add flaps!

end # module Dreadnought
