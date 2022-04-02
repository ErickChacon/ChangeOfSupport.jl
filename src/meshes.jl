# Methods for structs defined in Meshes.jl and additional meshes-related structs.

for filename in ["cartesiangrid.jl", "rectilineargrid.jl", "regularknots.jl"]
    include(joinpath("meshes", filename))
end

