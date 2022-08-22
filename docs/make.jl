using Documenter
using DynamicElectricSheath

DocMeta.setdocmeta!(DynamicElectricSheath, :DocTestSetup, :(using DynamicElectricSheath); recursive=true)

makedocs(;
    sitename = "DynamicElectricSheath.jl",
    repo="https://github.com/JuliaVlasov/DynamicElectricSheath.jl/blob/{commit}{path}#{line}",
    format = Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaVlasov.github.io/DynamicElectricSheath.jl",
        edit_link="main",
        assets=String[]),
    modules = [DynamicElectricSheath],
    pages = ["index.md", "example.md"]
)


deploydocs(;
    repo="github.com/JuliaVlasov/DynamicElectricSheath.jl",
    devbranch="main",
)
