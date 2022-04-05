# Methods for structs defined in Meshes.jl and additional meshes-related structs.

for filename in ["abstract.jl", "gmrf.jl", "rgmrf.jl", "cgmrf.jl"]
    include(joinpath("gmrf", filename))
end

