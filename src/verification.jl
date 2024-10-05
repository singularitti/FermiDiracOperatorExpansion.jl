using Roots: Newton, find_zero

export estimate_mu, rescale_mu

function estimate_mu(Nocc, 𝛆, β, μ₀; atol=1e-8, maxiters=50, verbose=true)
    g(μ) = Nocc - sum(fermi_dirac.(𝛆, μ, β))
    g′(μ) = sum(fermi_dirac_prime.(𝛆, μ, β))
    return find_zero((g, g′), μ₀, Newton(); atol=atol, maxiters=maxiters, verbose=verbose)
end

rescale_mu(μ′, α, μ) = (μ′ - oneunit(μ′) / 2) / α + μ
# rescale_mu(μ′, α, μ) = α * (μ′ - μ) + oneunit(μ′) / 2
