export fermi_dirac

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
