using FermiDiracOperatorExpansion
using BenchmarkTools
using GershgorinDiscs
using LinearAlgebra: Symmetric, diagm, diag, tr, eigvals

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

function estimate_mu(Nocc, 𝐇)
    Nocc = floor(Int, Nocc)
    diagonals = sort(diag(𝐇))
    HOMO, LUMO = diagonals[Nocc], diagonals[Nocc + 1]
    μ₀ = (HOMO + LUMO) / 2
    return estimate_mu(Nocc, diagonals, β, μ₀)
end
function compute_mu(Nocc, 𝐇)
    Nocc = floor(Int, Nocc)
    evals = eigvals(𝐇)
    HOMO, LUMO = evals[Nocc], evals[Nocc + 1]
    μ₀ = (HOMO + LUMO) / 2
    return estimate_mu(Nocc, evals, β, μ₀)
end

β = 40
μ = 0.6
order = 2^20
𝐇 = setup_hamiltonian(1000)

α = estimate_alpha(𝐇, μ)
α_exact = compute_alpha(𝐇, μ)

dm = density_matrix(𝐇, μ, α; order)
N = tr(dm)

# dm_α_exact = density_matrix(𝐇, μ, α_exact; order=order_α_exact)
# N_α_exact = tr(dm_α_exact)

dm_exact = fermi_dirac(𝐇, μ, β)
N_exact = tr(dm_exact)
