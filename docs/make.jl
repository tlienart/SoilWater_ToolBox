using Documenter
using SoilWater-ToolBox

makedocs(
    sitename = "SoilWater-ToolBox",
    format = Documenter.HTML(),
    modules = [SoilWater-ToolBox]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
