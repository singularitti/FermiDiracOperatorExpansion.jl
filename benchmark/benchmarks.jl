using FermiDiracOperatorExpansion
using BenchmarkTools
using GershgorinDiscs
using LinearAlgebra: Symmetric, diagm, diag, tr, eigvals
using Roots: Newton, find_zero

SUITE = BenchmarkGroup()
SUITE["rand"] = @benchmarkable rand(10)

# Write your benchmarks here.
function setup_hamiltonian(N, a=0.01)
    𝐇 = diagm(10.0 * (1:N))
    foreach(1:size(𝐇, 1)) do i
        foreach((i + 1):size(𝐇, 2)) do j
            𝐇[i, j] = exp(-a * (i - j)^2)  # Mimic a non-metallic system or a metallic system at ﬁnite temperature
        end
    end
    return Symmetric(𝐇)
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
    # return find_zero((g, g′), μ₀, Newton(); atol=1e-8, maxiters=50, verbose=true)
    return newton_raphson(g, g′, μ₀; tol=1e-10, maxiter=1000)
end
function compute_mu(𝐇, nocc)
    nocc = floor(Int, nocc)
    evals = eigvals(𝐇)
    HOMO, LUMO = evals[nocc], evals[nocc + 1]
    μ₀ = (HOMO + LUMO) / 2
    @show μ₀
    g(μ) = nocc - sum(fermi_dirac.(evals, μ, β))
    g′(μ) = sum(fermi_dirac_derivative.(evals, μ, β))
    # return find_zero((g, g′), μ₀, Newton(); atol=1e-8, maxiters=50, verbose=true)
    return newton_raphson(g, g′, μ₀; tol=1e-10, maxiter=1000)
end

function newton_raphson(
    f::Function,
    fprime::Function,
    x0::Number,
    args::Tuple=();
    tol::AbstractFloat=1e-8,
    maxiter::Integer=50,
    eps::AbstractFloat=1e-10,
)
    for _ in 1:maxiter
        yprime = fprime(x0, args...)
        if abs(yprime) < eps
            return x0
        end
        y = f(x0, args...)
        x1 = x0 - y / yprime
        if abs(x1 - x0) < tol
            @show x1 - x0
            return x1
        end
        x0 = x1
    end
    return error("Max iteration exceeded")
end

β = 1.1604441716111258
μ = 0.1
𝐇 = setup_hamiltonian(100)

# α = estimate_alpha(𝐇, μ)
# α_exact = compute_alpha(𝐇, μ)

# order = get_order(α, β)
# order_α_exact = get_order(α_exact, β)

# dm = density_matrix(𝐇, μ, α; order)
# N = tr(dm) / size(dm, 1)

# dm_α_exact = density_matrix(𝐇, μ, α_exact; order=order_α_exact)
# N_α_exact = tr(dm_α_exact) / size(dm, 1)

# dm_exact = fermi_dirac(𝐇, μ, β)
# N_exact = tr(dm_exact) / size(dm, 1)
