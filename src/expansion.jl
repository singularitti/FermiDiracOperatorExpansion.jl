export density_matrix, estimate_alpha, compute_alpha, normalize, expand

using GershgorinDiscs: eigvals_extrema
using LinearAlgebra: I, Diagonal, checksquare, eigen, eigvals
using OffsetArrays: OffsetVector, Origin

function expand(𝐗₀::AbstractMatrix, order=2048)
    𝐗₀ = collect(𝐗₀)
    checksquare(𝐗₀)  # See https://discourse.julialang.org/t/120556/2
    𝐗ᵢ = 𝐗₀  # i=0
    iterations = OffsetVector([𝐗ᵢ], Origin(0))
    foreach(1:ceil(log2(order))) do _  # Start from i+1
        𝐗ᵢ² = 𝐗ᵢ^2
        𝐀 = 2𝐗ᵢ² - 2𝐗ᵢ + I
        𝐗ᵢ = 𝐀 \ 𝐗ᵢ² # It is actually 𝐗ᵢ₊₁ = (𝐗ᵢ² + (I - 𝐗ᵢ)²)⁻¹ 𝐗ᵢ²
        push!(iterations, 𝐗ᵢ)
    end
    return iterations
end

function estimate_alpha(𝐇::AbstractMatrix, μ)
    λₘᵢₙ, λₘₐₓ = eigvals_extrema(𝐇)
    return minimum((inv(μ - λₘᵢₙ), inv(λₘₐₓ - μ))) / 2
end

function compute_alpha(𝐇::AbstractMatrix, μ)
    λₘᵢₙ, λₘₐₓ = extrema(eigvals(𝐇))
    return minimum((inv(μ - λₘᵢₙ), inv(λₘₐₓ - μ))) / 2
end

normalize(𝐇::AbstractMatrix, μ, α=estimate_alpha(𝐇, μ)) = α * (𝐇 - μ * I) + I / 2

function density_matrix(𝐇::AbstractMatrix, μ, α=estimate_alpha(𝐇, μ), order=2048)
    𝐗₀ = normalize(𝐇, μ, α)
    iterations = expand(𝐗₀, order)
    𝐗ₙ = last(iterations)
    return I - 𝐗ₙ
end
