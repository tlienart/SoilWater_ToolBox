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
	include("Tool.jl")
	include("Read.jl")
	include("Hydro\\HydroStruct.jl")
	include("Hydro\\HydroInitialize.jl")
	include("Hydro\\WaterRetentionCurve.jl")
	include("Hydro\\Kunsat.jl")
	include("Stats.jl")
	include("Hydro\\ObjectiveFunction_Hydro.jl")
	include("Hydro\\START_HydroParam.jl")
	include("Hydro\\HydroRelation.jl")
	if option.Infilt
		include("Infilt\\Diffusivity.jl")
		include("Infilt\\Sorptivity.jl")
		include("Infilt\\Best.jl")
		include("Infilt\\OptInfilt.jl")
		include("Infilt\\START_Infilt.jl")
	end

	include("Psd\\PsdStruct.jl")
	include("Psd\\PsdInitialize.jl")
	include("Psd\\PsdThetar.jl")
	include("Psd\\PsdFunc.jl")
	include("Psd\\PsdOpt.jl")
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

		if option.hydro.KunsatΨ
			K_KΨ, Ψ_KΨ, N_KΨ = read.KUNSATΨ(Id_Select, N_SoilSelect)
		else
			K_KΨ = zeros(Float64, N_SoilSelect,1)
			Ψ_KΨ = zeros(Float64, N_SoilSelect,1)
			N_KΨ = zeros(Float64, N_SoilSelect,1)
		end

		if option.Psd
			Rpart, ∑Psd, N_Psd, Φ_Psd = read.PSD(Id_Select, N_SoilSelect)
		else
			∑Psd = zeros(Float64, N_SoilSelect,1)
		end
		
		if option.Infilt
			Tinfilt, ∑Infilt, N_Infilt, infilt  = read.INFILTRATION(Id_Select, N_SoilSelect)
		end

		# Reinforcing the maximum of iSoil to simulate
		N_SoilSelect = min(N_SoilSelect, param.N_iSoil_Simulations)
	println("=== END  : READING === \n")


	if option.θΨ ≠ "No"
		println("=== START: DERIVING HYDRO PARAMETERS  ===")
			hydro = hydroParam.START_HYDROPARAM(N_SoilSelect, ∑Psd, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ)
		println("=== END  : DERIVING HYDRO PARAMETERS  === \n")
	else
		hydro = []
	end


	if option.Psd
		println("=== START: PSD MODEL  ===")
			psdparam, N_Psd, θ_Rpart, Ψ_Rpart, Psd = psd.START_PSD(N_SoilSelect, Ψ_θΨ, θ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, Rpart, ∑Psd, N_Psd, Φ_Psd, hydro)
		println("=== END  : PSD MODEL  === \n")
	else
		θ_Rpart = zeros(Float64, N_SoilSelect,1)
		Ψ_Rpart = zeros(Float64, N_SoilSelect,1)
	end

	
	if option.Infilt
		println("=== START: INFILTRATION  ===")
			∑Infilt_Best_HydroObs, ∑Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange, Tinfilt_Best_HydroObs = mainInfilt.START_INFILT(N_SoilSelect, Tinfilt, ∑Infilt, ∑Psd, N_Infilt, infilt, hydro)
		println("=== END  : INFILTRATION  === \n")
	end

	println("=== START: WRITING TABLE  ===")
		if option.θΨ ≠ "No"
			table.hydroParam.θΨK(Id_Select[1:N_SoilSelect], N_SoilSelect, hydro)
		end
		if option.Psd
			table.psd.PSD(Id_Select[1:N_SoilSelect], N_SoilSelect, psdparam)
		end  # if: name
	println("=== END  : WRITING TABLE  === \n")


	if option.Plot
		println("=== START: PLOTTING  ===")

		if option.θΨ ≠ "No" && option.hydro.Plot_θΨ
			plot.HYDROPARAM(Id_Select, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, N_SoilSelect, N_Psd, θ_Rpart, Ψ_Rpart, hydro)
		end # option.Plot_WaterRetentionCurve

		if option.Psd && option.psd.Plot_θr
			plot.PLOT_θr(∑Psd, N_SoilSelect, hydro, psdparam)
		end

		if option.Psd && option.psd.Plot_IMP_model
			plot.PLOT_IMP_model(Id_Select, Rpart, N_Psd, ∑Psd, Psd, N_SoilSelect, hydro, psdparam) 
		end
		# if option.Plot_BestLab && option.θΨ ≠ "No" && (option.infilt.OptimizeRun == "Run" ||  option.infilt.OptimizeRun == "RunOpt") && optiotestn.infiltration.Model=="Simplified" && option.infilt.SeIni_Range	

		# 	# plot.BEST_LAB_SEINIRANGE(Id_Select, ∑Infilt_Best_HydroObs_SeIniRange, N_SoilSelect, T_Best_HydroObs_SeIniRange)
		# end # option.Plot_BestLab

		# if option.Plot_BestLab && option.θΨ ≠ "No" && option.infilt.OptimizeRun  == "Run" || option.infilt.OptimizeRun  == "RunOpt" 
		# 	# plot.BEST_LAB(Id_Select, N_Infilt, N_SoilSelect, ∑Infilt_Best_HydroObs, Tinfilt_Best_HydroObs, Tinfilt, ∑Infilt)
		# end
	
		println("=== END: PLOTTING  === \n")
	end # if option.Plot
		
end  # function: START_TOOLBOX


println("\n\n===== START SOIL WATER TOOLBOX ==== \n")
	
		@time START_TOOLBOX()

println("==== END SOIL WATER TOOLBOX ===")