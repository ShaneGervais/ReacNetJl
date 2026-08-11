# Package extension implementing `_sparse_newton_solve` (see
# src/solver/sparse_solve.jl) with KLU.jl -- the sparse direct solver the
# feature spec recommends as primary (item #3) for this network structure.
#
# A fresh symbolic + numeric factorization is done on every call (one per
# Newton iteration); the sparsity *structure* never changes across a run
# (only nzval does, see `_cached_network_jacobian_sparse!`), so reusing the
# symbolic factorization across iterations via `KLU.klu!` is a real available
# speedup left for later -- this is the correct, simple version first.
module ReacNetJlKLUExt

using ReacNetJl
using KLU
using SparseArrays

import ReacNetJl: _sparse_newton_solve

function ReacNetJl._sparse_newton_solve(A::SparseMatrixCSC, b::Vector{Float64})
    F = klu(A)
    return F \ b
end

end # module
