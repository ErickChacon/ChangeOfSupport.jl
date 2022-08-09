# Methods for bsplines.

for filename in ["regular.jl", "basis.jl", "integral.jl", "nregular.jl"]
    include(joinpath("bsplines", filename))
end

