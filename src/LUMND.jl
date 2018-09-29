module LUMND
    
    include("../Metis/src/Metis.jl")

    export factorize

    using Printf

    include("util.jl")
    include("partition.jl")
    include("factorize.jl")
    include("solve.jl")

end
