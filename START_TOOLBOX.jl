##========================================================================================
##                                                                                      ##
##                                 Soil Water ToolBox                                   ##
##                                                                                      ##
##========================================================================================
using Suppressor

@suppress_err begin
	include("Option.jl")
	# Install packages to run program
		if option.DownloadPackage
			include("Packages.jl")
		end
	include("Path.jl")
	include("Cst.jl")
	include("Param.jl")
	include("Read.jl")
	include("HydroParam\\HydroStruct.jl")
	include("HydroParam\\WaterRetentionCurve.jl")
	include("Stats.jl")
	include("HydroParam\\Kunsat.jl")
	include("HydroParam\\ObjectiveFunction_Hydro.jl")
	include("HydroParam\\START_HydroParam.jl")
	include("Infilt\\Diffusivity.jl")
	include("Infilt\\Sorptivity.jl")
	include("HydroParam\\HydroRelation.jl")
	include("Infilt\\Best.jl")
	include("Infilt\\OptInfilt.jl")
	include("Infilt\\MAIN_Infilt.jl")
	include("Psd\\PsdThetar.jl")
	include("Psd\\PsdStruct.jl")
	include("Psd\\PsdInitiate.jl")
	include("Psd\\PsdFunc.jl")
	include("Psd\\START_PSD.jl")
	include("Table.jl")
	include("Plot.jl")
end


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
			Rpart, ∑Psd, N_Psd  = read.PSD(Id_Select, N_SoilSelect)
		else
			∑Psd = zeros(Float64, N_SoilSelect,1)
		end
		
		if option.Infiltration
			Tinfilt, ∑Infilt, N_Infilt, infilt  = read.INFILTRATION(Id_Select, N_SoilSelect)
		end

		# Reinforcing the maximum of iSoil to simulate
		N_SoilSelect = min(N_SoilSelect, param.N_iSoil_Simulations)
	println("=== END  : READING === \n")


	if option.θΨ ≠ "No"
		println("=== START: DERIVING HYDRO PARAMETERS  ===")
		Of, Of_θΨ, Of_Kunsat, hydro =  mainHydroParam.START_HYDROPARAM(N_SoilSelect, ∑Psd, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ)
		println("=== END  : DERIVING HYDRO PARAMETERS  === \n")
	else
		hydro = []
	end

	
	if option.Infiltration
		println("=== START: INFILTRATION  ===")
			∑Infilt_Best_HydroObs, ∑Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange, Tinfilt_Best_HydroObs = mainInfilt.MAIN_INFILT(N_SoilSelect, Tinfilt, ∑Infilt, ∑Psd, N_Infilt, infilt, hydro)
		println("=== END  : INFILTRATION  === \n")
	end


	if option.Psd
		println("=== START: PSD MODEL  ===")
			psdparam, θ_Rpart, Ψ_Rpart = psd.START_PSD(N_SoilSelect, Ψ_θΨ, θ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, Rpart, ∑Psd, N_Psd, hydro)
		println("=== END  : PSD MODEL  === \n")
	end


	println("=== START: TABLE  ===")
		if option.θΨ ≠ "No"
			table.θΨK(Id_Select[1:N_SoilSelect], Of, Of_θΨ, Of_Kunsat, N_SoilSelect, hydro)
		end
	println("=== END  : TABLE  === \n")


	if option.Plot
		println("=== START: PLOTTING  ===")
			if option.Plot_WaterRetentionCurve
				plot.HYDROPARAM(Id_Select, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, N_SoilSelect, hydro)
			end # option.Plot_WaterRetentionCurve

			if option.Plot_BestLab && option.θΨ ≠ "No" && (option.infiltration.OptimizeRun == "Run" ||  option.infiltration.OptimizeRun == "RunOpt") && option.infiltration.Model=="Simplified" && option.infiltration.SeIni_Range	

				plot.BEST_LAB_SEINIRANGE(Id_Select, ∑Infilt_Best_HydroObs_SeIniRange, N_SoilSelect, T_Best_HydroObs_SeIniRange)

			end # option.Plot_BestLab

			if option.Plot_BestLab && option.θΨ ≠ "No" && option.infiltration.OptimizeRun  == "Run" || option.infiltration.OptimizeRun  == "RunOpt" 
				plot.BEST_LAB(Id_Select, N_Infilt, N_SoilSelect, ∑Infilt_Best_HydroObs, Tinfilt_Best_HydroObs, Tinfilt, ∑Infilt)
			end

		println("=== END: PLOTTING  === \n")
		return
	end
		
end  # function: START_TOOLBOX


println("\n\n===== START SOIL WATER TOOLBOX ==== \n")
	
		@time START_TOOLBOX()

println("==== END SOIL WATER TOOLBOX ===")