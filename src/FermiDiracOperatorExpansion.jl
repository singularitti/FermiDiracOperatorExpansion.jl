module FermiDiracOperatorExpansion

export CG, NewtonSchulz
export expand

using ConjugateGradient: cg
using LinearAlgebra: I

abstract type Solver end
struct CG <: Solver end
struct NewtonSchulz <: Solver end

function expand(𝐗₀, ::CG; order=2048, cgiter=2000)
    M, N = size(𝐗₀)
    if M != N
        throw(DimensionMismatch("𝐗₀ must be a square matrix!"))
    end
    𝐗ᵢ = 𝐗₀  # i=0
    map(1:ceil(log2(order))) do _  # Start from i+1
        𝐗ᵢ² = 𝐗ᵢ^2
        𝐀 = 2𝐗ᵢ² - 2𝐗ᵢ + I
        𝐗ᵢ = splat(hcat)(
            map(1:N) do j
                𝐱 = 𝐗ᵢ[:, j]
                𝐛 = 𝐗ᵢ²[:, j]
                isconverged = false
                while isconverged
                    𝐱, _, isconverged = cg(𝐀, 𝐛, 𝐱; atol=eps(), maxiter=cgiter)
                end
                𝐱  # The jth column of 𝐗ᵢ₊₁
            end,
        )
        𝐗ᵢ
    end
    return nothing
end

end
