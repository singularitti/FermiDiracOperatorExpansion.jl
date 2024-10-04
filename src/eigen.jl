function homolumo_gap(𝐯, 𝐰, λₘᵢₙ, λₘₐₓ)
    map(zip(𝐯, 𝐰)) do v, w
        y₁ = (1 - sqrt(1 - 4v^2 / w)) / 2
        y₂ = (1 - sqrt(1 - 4v)) / 2
        y₃ = (1 + sqrt(1 - 4v)) / 2
        y₄ = (1 + sqrt(1 - 4v^2 / w)) / 2
        foreach(i:-1:1) do j
            if pj
                yₖ = sqrt(yₖ)
                yₖ = (yₖ - 1 + αⱼ) / αⱼ
            else
                yₖ = 1 - sqrt(1 - yₖ)
                yₖ = yₖ / αⱼ
            end
        end
        x₁ = min(x₁, y₁)
        x₂ = min(x₂, y₂)
        x₃ = max(x₃, y₃)
        x₄ = max(x₄, y₄)
    end
end
