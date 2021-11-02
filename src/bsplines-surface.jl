struct NRegularBsplines{N,T}
    lower::NTuple{N,T}
    upper::NTuple{N,T}
    df::NTuple{N, Int}
    order::Int
end

function basis(x, b::NRegularBsplines{N, T}) where {N, T}
    b = [RegularBsplines(b.lower[i], b.upper[i], b.df[i], b.order) for i in 1:N]
    bb = [basis(x[i], b[i]) for i in N:-1:1]
    N == 1 ? bb[1] : kron(bb...)
end

function basis(x::Array{T, N}, b::NRegularBsplines{N, T}) where {N, T}
    b = [RegularBsplines(b.lower[i], b.upper[i], b.df[i], b.order) for i in 1:N]
    bb = [basis(x[:, i], b[i]) for i in N:-1:1]
    N == 1 ? bb[1] : kron(bb...)
end

