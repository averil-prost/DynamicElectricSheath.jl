using Documenter
using FiniteDiffDynSheath

makedocs(
    sitename = "FiniteDiffDynSheath",
    format = Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true"),
    modules = [FiniteDiffDynSheath],
    pages = ["index.md", "example.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#



