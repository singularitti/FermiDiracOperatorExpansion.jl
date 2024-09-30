module FermiDiracOperatorExpansion

export CG, NewtonSchulz
export density_matrix, estimate_alpha, expand, fermi_dirac, get_order

using ConjugateGradient: cg
using GershgorinDiscs: eigvals_extrema
using LinearAlgebra: I, Diagonal, checksquare, eigen
using OffsetArrays: OffsetVector, Origin

# See https://github.com/JuliaMath/Roots.jl/blob/bf0da62/src/utils.jl#L9-L11
struct ConvergenceFailed
    msg::String
end

abstract type Solver end
Base.@kwdef struct CG <: Solver
    abstol::Float64 = 0
    maxiter::UInt64 = 2000
end
struct NewtonSchulz <: Solver end

function expand(𝐗₀::AbstractMatrix, solver::CG=CG(); order=2048)
    𝐗₀ = collect(𝐗₀)
    checksquare(𝐗₀)  # See https://discourse.julialang.org/t/120556/2
    𝐗ᵢ = 𝐗₀  # i=0
    iterations = OffsetVector([𝐗ᵢ], Origin(0))
    foreach(1:ceil(log2(order))) do _  # Start from i+1
        𝐗ᵢ² = 𝐗ᵢ^2
        𝐀 = 2𝐗ᵢ² - 2𝐗ᵢ + I
        𝐗ᵢ = splat(hcat)(
            map(zip(eachcol(𝐗ᵢ), eachcol(𝐗ᵢ²))) do (𝐱, 𝐛)
                𝐱, _, isconverged = cg(𝐀, 𝐛, 𝐱; atol=solver.atol, maxiter=solver.maxiter)
                if !isconverged
                    throw(ConvergenceFailed("CG did not converge! Increase `maxiter`!"))
                end
                𝐱  # Each column of 𝐗ᵢ₊₁
            end,
        )  # It is actually 𝐗ᵢ₊₁
        push!(iterations, 𝐗ᵢ)
    end
    return iterations
end

function estimate_alpha(𝐇::AbstractMatrix, μ)
    λₘᵢₙ, λₘₐₓ = eigvals_extrema(𝐇)
    return minimum((inv(μ - λₘᵢₙ), inv(λₘₐₓ - μ))) / 2
end

normalize(𝐇::AbstractMatrix, μ, α=estimate_alpha(𝐇, μ)) = α * (𝐇 - μ * I) + I / 2

function density_matrix(
    𝐇::AbstractMatrix, μ, α=estimate_alpha(𝐇, μ); solver::Solver=CG(), order=2048
)
    𝐗₀ = normalize(𝐇, μ, α)
    iterations = expand(𝐗₀, solver; order)
    𝐗ₙ = last(iterations)
    return I - 𝐗ₙ
end

get_order(α, β) = β / 4α

function fermi_dirac(ε, μ, β)
    η = exp((ε - μ) * β)
    return inv(oneunit(η) + η)
end
function fermi_dirac(𝐇::AbstractMatrix, μ, β)
    E = eigen(𝐇)
    Λ, V = E.values, E.vectors
    FD = fermi_dirac.(Λ, μ, β)
    return V * Diagonal(FD) * V'
end

end
