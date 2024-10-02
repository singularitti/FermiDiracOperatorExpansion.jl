using LinearAlgebra: norm, tr

export trace

function trace(X::AbstractMatrix)
    Δ = X - X^2
    v = norm(Δ, 2)  # Frobenius norm
    w = tr(Δ)
    return v, w
end
function trace(iterations::AbstractVector{<:AbstractMatrix{T}}) where {T}
    𝐯, 𝐰 = similar(iterations, T), similar(iterations, T)
    foreach(eachindex(iterations)) do i
        𝐯[i], 𝐰[i] = trace(iterations[i])
    end
    return 𝐯, 𝐰
end
