using Documenter
using FiniteDiffDynSheath

makedocs(
    sitename = "FiniteDiffDynSheath",
    format = Documenter.HTML(),
    modules = [FiniteDiffDynSheath]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
