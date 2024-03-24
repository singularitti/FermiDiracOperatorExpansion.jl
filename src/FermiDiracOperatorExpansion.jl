module FermiDiracOperatorExpansion

export CG, NewtonSchulz
export expand

using ConjugateGradient: cg
using LinearAlgebra: I
using OffsetArrays: OffsetVector, Origin

# See https://github.com/JuliaMath/Roots.jl/blob/bf0da62/src/utils.jl#L9-L11
struct ConvergenceFailed
    msg::String
end

abstract type Solver end
Base.@kwdef struct CG <: Solver
    atol::Float64 = eps()
    maxiter::UInt64 = 2000
end
struct NewtonSchulz <: Solver end

function expand(𝐗₀::AbstractMatrix, ::CG; order=2048, cgiter=2000)
    M, N = size(𝐗₀)
    if M != N
        throw(DimensionMismatch("𝐗₀ must be a square matrix!"))
    end
    𝐗ᵢ = 𝐗₀  # i=0
    iterations = OffsetVector([𝐗ᵢ], Origin(0))
    foreach(1:ceil(log2(order))) do _  # Start from i+1
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
        push!(iterations, 𝐗ᵢ)
    end
    return iterations
end

end
