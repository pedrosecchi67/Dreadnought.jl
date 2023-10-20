using Dreadnought

const time_benchmarks = true

include("n0012.jl")

if time_benchmarks
    @info "Running time benchmarks"

    include("time_benchmark/time_benchmark.jl")
end

include("23012_airfoil/n23012.jl")
include("2412_airfoil/n2412.jl")
include("4412_airfoil/n4412.jl")
include("6712_airfoil/n6712.jl")
include("8318_airfoil/n8318.jl")
include("nlf10215_airfoil/nlf10215.jl")
include("4415_s20deg_wing/n4415.jl")
