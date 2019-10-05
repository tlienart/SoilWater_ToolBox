include("Option.jl")
# Install packages to run program
	if option.DownloadPackage
		include("Packages.jl")
	end
include("Path.jl")
include("Read.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : START_TOOLBOX
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function START_TOOLBOX()

	println("=== START READING === \n")
		# Selecting soils of interest 
		Id_Select, Id_True, N_SoilSelect = read.ID()

		if option.θΨ ≠ "No"
			θ_θΨ, Ψ_θΨ, N_θΨ = read.θΨ(Id_Select, N_SoilSelect)
		end

		if option.KunsatΨ
			K_KΨ, Ψ_KΨ, N_KΨ = read.KUNSATΨ(Id_Select, N_SoilSelect)
		end

		if option.Psd
			Diameter, ∑Psd, N_Psd  = read.PSD(Id_Select, N_SoilSelect)
		end
		
		if option.Infiltration
			T, ∑Inf, N_Inf  = read.INFILTRATION(Id_Select, N_SoilSelect)
		end
	println("=== END READING ===")
		
end  # function: START_TOOLBOX


println("===== START SOIL WATER TOOLBOX ====")
	@time START_TOOLBOX()
println("==== END SOIL WATER TOOLBOX ===")