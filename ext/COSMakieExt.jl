module COSMakieExt

using Meshes
using ChangeOfSupport
using Makie

# Extension of rangebars for CartesianGrid{1}

function Makie.convert_arguments(P::Type{<:Makie.Rangebars}, d::CartesianGrid{1}, x::AbstractVector)
    marks = range(d)[1]
    Makie.convert_arguments(P, x, marks[1:end-1], marks[2:end])
end

end
