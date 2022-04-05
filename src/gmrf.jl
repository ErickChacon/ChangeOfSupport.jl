# Methods for abstract gmrf struct and subtypes as gmrf, rgmrf, cgmrf.

for filename in ["abstract.jl", "gmrf.jl", "rgmrf.jl", "cgmrf.jl"]
    include(joinpath("gmrf", filename))
end

