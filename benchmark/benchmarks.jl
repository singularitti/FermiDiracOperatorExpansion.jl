using FermiDiracOperatorExpansion
using BenchmarkTools
using LinearAlgebra: Diagonal, Symmetric

SUITE = BenchmarkGroup()
SUITE["rand"] = @benchmarkable rand(10)

# Write your benchmarks here.

const kBT = 0.25
const μ = 0.1

function setup_hamiltonian(N, a=0.01)
    𝐇 = Diagonal(10 * rand(N))
    foreach(1:N) do i
        foreach(1:N) do j
            𝐇[i, j] = exp(-a * (i - j)^2)  # Mimic a non-metallic system or a metallic system at ﬁnite temperature
        end
    end
    return Symmetric(𝐇)
end
