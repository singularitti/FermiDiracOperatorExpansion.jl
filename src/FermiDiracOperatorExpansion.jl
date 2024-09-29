module FermiDiracOperatorExpansion

export CG, NewtonSchulz
export densitymatrix, estimate_alpha, expand

using ConjugateGradient: cg
using GershgorinDiscs: eigvals_extrema
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

function expand(𝐗₀::AbstractMatrix, solver::CG; order=2048)
    𝐗₀ = collect(𝐗₀)
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
                𝐱, _, isconverged = cg(𝐀, 𝐛, 𝐱; atol=solver.atol, maxiter=solver.maxiter)
                if !isconverged
                    throw(ConvergenceFailed("CG did not converge! Increase `maxiter`!"))
                end
                𝐱  # The jth column of 𝐗ᵢ₊₁
            end,
        )
        push!(iterations, 𝐗ᵢ)
    end
    return iterations
end

function estimate_alpha(𝐇::AbstractMatrix, mu=1 / 2)
    λₘᵢₙ, λₘₐₓ = eigvals_extrema(𝐇)
    return minimum((inv(mu - λₘᵢₙ), inv(λₘₐₓ - mu))) / 2
end

normalize(𝐇::AbstractMatrix; mu=1 / 2, alpha=estimate_alpha(𝐇, mu)) =
    alpha * (𝐇 - mu * I) + I / 2

function densitymatrix(
    𝐇::AbstractMatrix, solver::Solver; mu=1 / 2, alpha=estimate_alpha(𝐇, mu), order=2048
)
    iterations = expand(normalize(𝐇; mu, alpha), solver; order)
    𝐗ₙ = last(iterations)
    return I - 𝐗ₙ
end

end
