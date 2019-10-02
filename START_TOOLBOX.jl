include("Option.jl")
# Install packages to run program
	if option.Packages
		include("Packages.jl")
	end
include("Path.jl")

