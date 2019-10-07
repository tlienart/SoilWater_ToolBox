include("Option.jl")
# Install packages to run program
	if option.DownloadPackage
		include("Packages.jl")
	end
include("Path.jl")
include("Param.jl")
include("Read.jl")
include("Psd\\PsdThetar.jl")
include("HydroParam\\WaterRetentionCurve.jl")
include("Stats.jl")
include("HydroParam\\Kunsat.jl")
include("HydroParam\\ObjectiveFunction_Hydro.jl")
include("HydroParam\\MAINhydroParam.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : START_TOOLBOX
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function START_TOOLBOX()

	println("=== START: READING ===")
		# Selecting soils of interest 
		Id_Select, Id_True, N_SoilSelect = read.ID()

		if option.θΨ ≠ "No"
			θ_θΨ, Ψ_θΨ, N_θΨ = read.θΨ(Id_Select, N_SoilSelect)
		end

		if option.KunsatΨ
			K_KΨ, Ψ_KΨ, N_KΨ = read.KUNSATΨ(Id_Select, N_SoilSelect)
		else
			K_KΨ = zeros(Float64, N_SoilSelect,1)
			Ψ_KΨ = zeros(Float64, N_SoilSelect,1)
			N_KΨ = zeros(Float64, N_SoilSelect,1)
		end

		if option.Psd
			Diameter, ∑Psd, N_Psd  = read.PSD(Id_Select, N_SoilSelect)
		else
			∑Psd = zeros(Float64, N_SoilSelect,1)
		end
		
		if option.Infiltration
			T, ∑Inf, N_Inf  = read.INFILTRATION(Id_Select, N_SoilSelect)
		end
	println("=== END  : READING === \n")


	if option.θΨ ≠ "No"
		println("=== START: DERIVING HYDRO PARAMETERS  ===")

			hydro =  mainHydroParam.MAIN_HYDROPARAM(N_SoilSelect, ∑Psd, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ)
		
		println("=== END  : DERIVING HYDRO PARAMETERS  === \n")
	end



		
end  # function: START_TOOLBOX


println("===== START SOIL WATER TOOLBOX ==== \n")
	@time START_TOOLBOX()
println("==== END SOIL WATER TOOLBOX ===")