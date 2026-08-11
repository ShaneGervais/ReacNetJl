# _sparse_newton_solve hook, implemented by the KLU package extension.

# Solve `A * x = b` for the backward-Euler Newton correction when
# `jacobian=:sparse` (see backward_euler.jl), where `A = I - dt*J` is built
# from the network's precomputed sparsity pattern (feature spec Tier 0 #2/#3).
# This function is provided by a package extension: install and load the
# solver package first with `import Pkg; Pkg.add("KLU")` and `using KLU`.
#
# Deliberately untyped (unlike the extension's `SparseMatrixCSC`/`Vector{Float64}`
# method): Julia extensions *add* methods to a stub, they cannot narrow/replace
# one with an identical signature -- precompiling the extension errors as
# "method overwriting" if this stub's signature exactly matches the real
# implementation's. Keeping the stub generic so the extension's method is
# strictly more specific (and therefore always preferred once loaded) is the
# same pattern `solve_network_fbdf`'s stub in fbdf.jl uses.
function _sparse_newton_solve(A, b)
    error("jacobian=:sparse requires the KLU extension; run `import Pkg; Pkg.add(\"KLU\")` and add `using KLU` before calling it")
end
