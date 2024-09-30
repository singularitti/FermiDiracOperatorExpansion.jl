using FermiDiracOperatorExpansion
using BenchmarkTools
using LinearAlgebra: Symmetric, diagm, tr

SUITE = BenchmarkGroup()
SUITE["rand"] = @benchmarkable rand(10)

# Write your benchmarks here.
function setup_hamiltonian(N, a=0.01)
    𝐇 = diagm(10 * rand(N))
    foreach(1:size(𝐇, 1)) do i
        foreach(1:size(𝐇, 2)) do j
            𝐇[i, j] = exp(-a * (i - j)^2)  # Mimic a non-metallic system or a metallic system at ﬁnite temperature
        end
    end
    return Symmetric(𝐇)
end

β = 4
μ = 0.1
𝐇 = setup_hamiltonian(100)
α = estimate_alpha(𝐇, μ)
order = get_order(α, β)
dm = density_matrix(𝐇, μ, α; order)
N = tr(dm)
dm_exact = fermi_dirac(𝐇, μ, β)
N_exact = tr(dm_exact)
