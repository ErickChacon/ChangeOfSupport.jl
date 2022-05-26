
# check xtx
sum(abs.(out - out2))

# check cholesky
using LinearAlgebra
ch1 = cholesky(out)
ch2 = cholesky(out2)

ch1.L
ch2.L |> sparse |> Matrix
ch2.p
# cholesky(sparse(A), perm=1:3)
P = sparse(1:50, ch2.p, ones(50))
L = sparse(ch2.L)
sum(abs.(P'*L*L'*P - out2))
sum(abs.(L*L' - out2[ch2.p, ch2.p]))

sum(abs.(L' \ collect(1:50) - ch2.U \ collect(1:50)))
sum(abs.((L'*P) \ collect(1:50) - ch2.UP \ collect(1:50)))
sum(abs.(ch2.U \ collect(1:50) - ch2.UP \ collect(1:50)))

sparse(out2.L) |> Matrix

histogram(ch1.U \ ones(50), nbins = 1000)
histogram(ch2.U \ ones(50), nbins = 1000)

histogram(ch2.UP \ collect(1:50), nbins = 1000)
histogram(ch2.U \ collect(1:50), nbins = 1000)
sum(abs.(sort(ch2.UP \ collect(1:50)) - sort(ch2.U \ collect(1:50))))

histogram(out.L \ collect(1:50), nbins = 20)
histogram(out2.PtL \ collect(1:50), nbins = 20)

