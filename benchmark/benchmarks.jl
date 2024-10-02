using FermiDiracOperatorExpansion
using BenchmarkTools
using GershgorinDiscs
using LinearAlgebra: Symmetric, diagm, diag, tr, eigvals
using Roots: Newton, find_zero

SUITE = BenchmarkGroup()
SUITE["rand"] = @benchmarkable rand(10)

# Write your benchmarks here.
function setup_hamiltonian(N, a=0.01)
    𝐇 = diagm(10.0 * rand(N))
    foreach(1:size(𝐇, 1)) do i
        foreach((i + 1):size(𝐇, 2)) do j
            𝐇[i, j] = exp(-a * (i - j)^2)  # Mimic a non-metallic system or a metallic system at ﬁnite temperature
        end
    end
    return Symmetric(𝐇)
end
function setup_hamiltonian2(N)
    𝐇 = zeros(N, N)
    foreach(1:size(𝐇, 1)) do i
        foreach((i + 1):size(𝐇, 2)) do j
            𝐇[i, j] = exp(-0.0005abs(i - j) / 2) * sin(i + j)
        end
    end
    return Symmetric(𝐇)
end
function setup_hamiltonian3(N)
    return 100 * diagm(sort(rand(N)))
end

function fermi_dirac_derivative(ε, μ, β)
    fd = fermi_dirac(ε, μ, β)
    return -β * fd * (oneunit(fd) - fd)
end

function estimate_mu(𝐇, nocc)
    nocc = floor(Int, nocc)
    diagonal = sort(diag(𝐇))
    HOMO, LUMO = diagonal[nocc], diagonal[nocc + 1]
    μ₀ = (HOMO + LUMO) / 2
    @show μ₀
    g(μ) = nocc - sum(fermi_dirac.(diagonal, μ, β))
    g′(μ) = sum(fermi_dirac_derivative.(diagonal, μ, β))
    return find_zero((g, g′), μ₀, Newton(); atol=1e-8, maxiters=50, verbose=true)
end
function compute_mu(𝐇, nocc)
    nocc = floor(Int, nocc)
    evals = eigvals(𝐇)
    HOMO, LUMO = evals[nocc], evals[nocc + 1]
    μ₀ = (HOMO + LUMO) / 2
    @show μ₀
    g(μ) = nocc - sum(fermi_dirac.(evals, μ, β))
    g′(μ) = sum(fermi_dirac_derivative.(evals, μ, β))
    return find_zero((g, g′), μ₀, Newton(); atol=1e-8, maxiters=50, verbose=true)
end

β = 1.1604441716111258
μ = 0.1
𝐇 = setup_hamiltonian(1000)

α = estimate_alpha(𝐇, μ)
α_exact = compute_alpha(𝐇, μ)

dm = density_matrix(𝐇, μ, α; order)
N = tr(dm)

# dm_α_exact = density_matrix(𝐇, μ, α_exact; order=order_α_exact)
# N_α_exact = tr(dm_α_exact)

dm_exact = fermi_dirac(𝐇, μ, β)
N_exact = tr(dm_exact)
