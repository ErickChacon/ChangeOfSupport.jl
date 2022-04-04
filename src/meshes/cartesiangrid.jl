## Basic methods for CartesianGrid

function Base.range(x::CartesianGrid{Dim}) where {Dim}
    lower = coordinates(minimum(x))
    upper = coordinates(maximum(x))
    step = spacing(x)
    [(lower[i]):(step[i]):(upper[i]) for i in 1:Dim]
end

function centroids(x::CartesianGrid{Dim}) where {Dim}
    gridknots = range(x)
    [knots[1:end-1] .+ 0.5 * step(knots) for knots in gridknots]
end

function centroidsmat(x::CartesianGrid{Dim}) where {Dim}
    Matrix(transpose(reduce(hcat, coordinates.(centroid.(x)))))
end

## Adjacency sparse matrix of CartesianGrid (temporal, spatial)

function adjacency(g::CartesianGrid{1}; order = 1, cyclic = false)
    n = nelements(g)

    # neighbors to the right (→)
    if cyclic
        Ir = 1:n
        Jr = mod1.(Ir .+ order, n)
    else
        Ir = 1:(n - order)
        Jr = (1 + order):n
    end

    # sparse adjacency matrix
    A = sparse(Ir, Jr, true, n, n)
    A .| A'
end

function adjacency(g::CartesianGrid{2}; order = 1, cyclic = false)
    n1, n2 = size(g)
    n = nelements(g)

    # function to convert i,j cel to k-index (CartesianGrid)
    function ij_to_k(i, j, n1, n2)
        (j - 1) * n1 + i
    end

    if cyclic
        # neighbors to the right (→)
        Ir = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:n2]
        Jr = [ij_to_k(mod1(i + order, n1), j, n1, n2) for i = 1:n1 for j = 1:n2]
        # neighbors to the top (↑)
        It = Ir
        Jt = [ij_to_k(i, mod1(j + order, n2), n1, n2) for i = 1:n1 for j = 1:n2]
    else
        # neighbors to the right (→)
        Ir = [ij_to_k(i, j, n1, n2) for i = 1:(n1-order) for j = 1:n2]
        Jr = [ij_to_k(i + order, j, n1, n2) for i = 1:(n1-order) for j = 1:n2]
        # neighbors to the top (↑)
        It = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:(n2-order)]
        Jt = [ij_to_k(i, j + order, n1, n2) for i = 1:n1 for j = 1:(n2-order)]
    end

    # sparse adjacency matrix
    A = sparse(vcat(Ir, It), vcat(Jr, Jt), true, n, n)
    A .| A'
end

## Differences and structure for GMRF

function difference(g::CartesianGrid{1}; order = 1, cyclic = false)
    n = nelements(g)

    if order == 1 && !cyclic
        return spdiagm(n-1, n, 0 => fill(-1, n-1), 1 => fill(1, n-1))
    end

    if order == 1 && cyclic
        return spdiagm(0 => fill(-1, n), 1 => fill(1, n-1), -(n-1) => fill(1, 1))
    end

    if order == 2 && !cyclic
        return spdiagm(n-2, n,
                       0 => fill(1, n-2), 1 => fill(-2, n-2), 2 => fill(1, n-2))
    end

    if order == 2 && cyclic
        return spdiagm(0 => fill(1, n), 1 => fill(-2, n-1), 2 => fill(1, n-2),
                       -(n-2) => fill(1, 2), -(n-1) => fill(-2, 1))
    end

    throw(ErrorException("not implemented"))
end

function difference(g::CartesianGrid{2}; order = 1, cyclic = false)
    n1, n2 = size(g)
    n = nelements(g)

    # function to convert i,j cel to k-index (CartesianGrid)
    function ij_to_k(i, j, n1, n2)
        (j - 1) * n1 + i
    end

    if order == 1
        if cyclic
            # neighbors to the right (→)
            Ir = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:n2]
            Jr = [ij_to_k(mod1(i + 1, n1), j, n1, n2) for i = 1:n1 for j = 1:n2]
            # neighbors to the top (↑)
            It = Ir
            Jt = [ij_to_k(i, mod1(j + 1, n2), n1, n2) for i = 1:n1 for j = 1:n2]
        else
            # neighbors to the right (→)
            Ir = [ij_to_k(i, j, n1, n2) for i = 1:(n1-1) for j = 1:n2]
            Jr = [ij_to_k(i + 1, j, n1, n2) for i = 1:(n1-1) for j = 1:n2]
            # neighbors to the top (↑)
            It = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:(n2-1)]
            Jt = [ij_to_k(i, j + 1, n1, n2) for i = 1:n1 for j = 1:(n2-1)]
        end
        return sparse(
            repeat(1:(length(Ir) + length(It)), 2),
            vcat(Ir, It, Jr, Jt),
            repeat([1, -1], inner = length(Ir) + length(It))
           )
    end

    if order == 2 && cyclic
        Ir = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:n2]
        # neighbors to the right (→)
        Jr = [ij_to_k(mod1(i + 1, n1), j, n1, n2) for i = 1:n1 for j = 1:n2]
        # neighbors to the top (↑)
        Jt = [ij_to_k(i, mod1(j + 1, n2), n1, n2) for i = 1:n1 for j = 1:n2]
        # neighbors to the left (←)
        Jl = [ij_to_k(mod1(i - 1, n1), j, n1, n2) for i = 1:n1 for j = 1:n2]
        # neighbors to the bottom (↓)
        Jb = [ij_to_k(i, mod1(j - 1, n2), n1, n2) for i = 1:n1 for j = 1:n2]
        return sparse(
            repeat(1:length(Ir), 5),
            vcat(Ir, Jr, Jt, Jl, Jb),
            repeat([-4, 1, 1, 1, 1], inner = length(Ir))
           )
    end

    if order == 2 && !cyclic
        Iir = [ij_to_k(i, j, n1, n2) for i = 2:(n1-1) for j = 2:(n2-1)]
        # neighbors to the right (→)
        Jir = [ij_to_k(i + 1, j, n1, n2) for i = 2:(n1-1) for j = 2:(n2-1)]
        # neighbors to the top (↑)
        Jit = [ij_to_k(i, j + 1, n1, n2) for i = 2:(n1-1) for j = 2:(n2-1)]
        # neighbors to the left (←)
        Jil = [ij_to_k(i - 1, j, n1, n2) for i = 2:(n1-1) for j = 2:(n2-1)]
        # neighbors to the bottom (↓)
        Jib = [ij_to_k(i, j - 1, n1, n2) for i = 2:(n1-1) for j = 2:(n2-1)]
        # sides
        Isl = [ij_to_k(1, j, n1, n2) for j = 2:(n2-1)]
        Jsl1 = [ij_to_k(1 + 1, j, n1, n2) for j = 2:(n2-1)]
        Jsl2 = [ij_to_k(1, j + 1, n1, n2) for j = 2:(n2-1)]
        Jsl3 = [ij_to_k(1, j - 1, n1, n2) for j = 2:(n2-1)]
        Isr = [ij_to_k(n1, j, n1, n2) for j = 2:(n2-1)]
        Jsr1 = [ij_to_k(n1 - 1, j, n1, n2) for j = 2:(n2-1)]
        Jsr2 = [ij_to_k(n1, j + 1, n1, n2) for j = 2:(n2-1)]
        Jsr3 = [ij_to_k(n1, j - 1, n1, n2) for j = 2:(n2-1)]
        Ist = [ij_to_k(i, n2, n1, n2) for i = 2:(n1-1)]
        Jst1 = [ij_to_k(i + 1, n2, n1, n2) for i = 2:(n1-1)]
        Jst2 = [ij_to_k(i - 1, n2, n1, n2) for i = 2:(n1-1)]
        Jst3 = [ij_to_k(i, n2 - 1, n1, n2) for i = 2:(n1-1)]
        Isb = [ij_to_k(i, 1, n1, n2) for i = 2:(n1-1)]
        Jsb1 = [ij_to_k(i + 1, 1, n1, n2) for i = 2:(n1-1)]
        Jsb2 = [ij_to_k(i - 1, 1, n1, n2) for i = 2:(n1-1)]
        Jsb3 = [ij_to_k(i, 1 + 1, n1, n2) for i = 2:(n1-1)]
        Is = vcat(Isl, Isr, Ist, Isb)
        Js = vcat(Jsl1, Jsr1, Jst1, Jsb1,
                  Jsl2, Jsr2, Jst2, Jsb2,
                  Jsl3, Jsr3, Jst3, Jsb3)
        # corners
        Ic1 = ij_to_k(1, 1, n1, n2)
        Jc1 = [ij_to_k(1 + 1, 1, n1, n2), ij_to_k(1, 1 + 1, n1, n2)]
        Ic2 = ij_to_k(1, n2, n1, n2)
        Jc2 = [ij_to_k(1 + 1, n2, n1, n2), ij_to_k(1, n2 - 1, n1, n2)]
        Ic3 = ij_to_k(n1, 1, n1, n2)
        Jc3 = [ij_to_k(n1 - 1, 1, n1, n2), ij_to_k(n1, 1 + 1, n1, n2)]
        Ic4 = ij_to_k(n1, n2, n1, n2)
        Jc4 = [ij_to_k(n1 - 1, n2, n1, n2), ij_to_k(n1, n2 - 1, n1, n2)]
        Ic = vcat(Ic1, Ic2, Ic3, Ic4)
        Jc = vcat(Jc1[1], Jc2[1], Jc3[1], Jc4[1], Jc1[2], Jc2[2], Jc3[2], Jc4[2])
        return sparse(
            vcat(repeat(1:length(Iir), 5),
                 repeat((1:length(Is)) .+ length(Iir), 4),
                 repeat((1:length(Ic)) .+ length(Iir) .+ length(Is), 3),
                ),
            vcat(Iir, Jir, Jit, Jil, Jib, Is, Js, Ic, Jc),
            vcat(repeat([-4, 1, 1, 1, 1], inner = length(Iir)),
                 repeat([-3, 1, 1, 1], inner = length(Is)),
                 repeat([-2, 1, 1], inner = length(Ic)),
                )
           )
    end

    throw(ErrorException("not implemented"))
end

function structure(g::CartesianGrid; δ = 0, order = 1, cyclic = false)
    D = difference(g, order = order, cyclic = cyclic)
    D'D + δ * I
end

