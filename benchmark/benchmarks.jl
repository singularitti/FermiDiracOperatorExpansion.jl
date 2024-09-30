using FermiDiracOperatorExpansion
using BenchmarkTools
using GershgorinDiscs
using LinearAlgebra: Symmetric, diagm, tr

SUITE = BenchmarkGroup()
SUITE["rand"] = @benchmarkable rand(10)

# Write your benchmarks here.
function setup_hamiltonian(N, a=0.01)
    𝐇 = diagm(10 * rand(N))
    foreach(1:size(𝐇, 1)) do i
        foreach((i + 1):size(𝐇, 2)) do j
            𝐇[i, j] = exp(-a * (i - j)^2)  # Mimic a non-metallic system or a metallic system at ﬁnite temperature
        end
    end
    return Symmetric(𝐇)
end

function estimate_mu(β, N, λₘᵢₙ, λₘₐₓ)
    η = exp(β * N)
    return inv(β) * log((η - 1) / (exp(-β * λₘᵢₙ) - η * exp(-β * λₘₐₓ)))
end

β = 4
μ = 0.1
𝐇 = setup_hamiltonian(1000)

α = estimate_alpha(𝐇, μ)
α_exact = compute_alpha(𝐇, μ)

order = get_order(α, β)
order_α_exact = get_order(α_exact, β)

dm = density_matrix(𝐇, μ, α; order)
N = tr(dm)
mu = estimate_mu(β, N, eigvals_extrema(𝐇)...)

dm_α_exact = density_matrix(𝐇, μ, α_exact; order=order_α_exact)
N_α_exact = tr(dm_α_exact)

dm_exact = fermi_dirac(𝐇, μ, β)
N_exact = tr(dm_exact)
mu_exact = estimate_mu(β, N_exact, extrema(eigvals(𝐇))...)
