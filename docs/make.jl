using FermiDiracOperatorExpansion
using Documenter

DocMeta.setdocmeta!(FermiDiracOperatorExpansion, :DocTestSetup, :(using FermiDiracOperatorExpansion); recursive=true)

makedocs(;
    modules=[FermiDiracOperatorExpansion],
    authors="singularitti <singularitti@outlook.com> and contributors",
    sitename="FermiDiracOperatorExpansion.jl",
    format=Documenter.HTML(;
        canonical="https://singularitti.github.io/FermiDiracOperatorExpansion.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/singularitti/FermiDiracOperatorExpansion.jl",
    devbranch="main",
)
