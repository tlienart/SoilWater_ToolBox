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
			include("Infilt\\Sorptivity.jl")
			include("Infilt\\Best.jl")
			include("Infilt\\OptInfilt.jl")
			include("Infilt\\START_Infilt.jl")
		end
		if option.Psd
			include("Psd\\PsdStruct.jl")
			include("Psd\\PsdInitialize.jl")
		end
	include("Psd\\PsdThetar.jl")
		if option.Psd
			include("Psd\\PsdFunc.jl")
			include("Psd\\PsdOpt.jl")
			include("Psd\\START_PSD.jl")
		end
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
		end

		if option.Psd
			Rpart, ∑Psd, N_Psd, Φ_Psd = read.PSD(Id_Select, N_SoilSelect)
		else
			∑Psd = zeros(Float64, N_SoilSelect,1)
		end
		
		if option.Infilt
			Tinfilt, ∑Infilt, N_Infilt, infiltParam  = read.INFILTRATION(Id_Select, N_SoilSelect)
		end

		# Reinforcing the maximum of iSoil to simulate
		N_SoilSelect = min(N_SoilSelect, param.N_iSoil_Simulations)
	println("=== END  : READING === \n")


	if option.θΨ ≠ "No"
		println("=== START: DERIVING HYDRO PARAMETERS  ===")
			# INITIALIZES HYDRAULIC PARAMETERS STRUCT INDEPENDENTLY OF THE SELECTED MODEL
			hydro = hydroStruct.HYDROSTRUCT(N_SoilSelect)

			if option.hydro.KunsatΨ
				# Structure of hydro
				hydro = hydroParam.START_HYDROPARAM(N_SoilSelect=N_SoilSelect, ∑Psd=∑Psd,  θ_θΨ=θ_θΨ, Ψ_θΨ=Ψ_θΨ, N_θΨ=N_θΨ, K_KΨ=K_KΨ, Ψ_KΨ=Ψ_KΨ, N_KΨ=N_KΨ, hydro=hydro, optionHydro=option.hydro)
			else
				hydro = hydroParam.START_HYDROPARAM(N_SoilSelect=N_SoilSelect, ∑Psd=∑Psd,  θ_θΨ=θ_θΨ, Ψ_θΨ=Ψ_θΨ, N_θΨ=N_θΨ, hydro=hydro, optionHydro=option.hydro)
			end
		println("=== END  : DERIVING HYDRO PARAMETERS  === \n")
	else
		hydro = []
	end


	if option.Psd
		println("=== START: PSD MODEL  ===")
			# Structure of psdHydro
			psdHydro = hydroStruct.HYDROSTRUCT(N_SoilSelect)

			psdParam, N_Psd, θ_Rpart, Ψ_Rpart, Psd, psdHydro = psd.START_PSD(N_SoilSelect, Ψ_θΨ, θ_θΨ, N_θΨ, Rpart, ∑Psd, N_Psd, Φ_Psd, hydro, psdHydro)

		if  option.psd.HydroParam
			psdHydro = hydroParam.START_HYDROPARAM(N_SoilSelect=N_SoilSelect, ∑Psd=∑Psd, θ_θΨ=θ_Rpart, Ψ_θΨ=Ψ_Rpart, N_θΨ=N_Psd, hydro=psdHydro, optionHydro=option.psd)
		end

		println("=== END  : PSD MODEL  === \n")
	else
		θ_Rpart = zeros(Float64, N_SoilSelect,1)
		Ψ_Rpart = zeros(Float64, N_SoilSelect,1)
	end

	
	if option.Infilt
		println("=== START: INFILTRATION  ===")
		# ∑Infilt_Best_HydroObs, ∑Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange, Tinfilt_Best_HydroObs =
		 infilt.START_INFILTRATION(N_SoilSelect, Tinfilt, ∑Infilt, ∑Psd, N_Infilt, infiltParam, hydro)
		println("=== END  : INFILTRATION  === \n")
	end

	println("=== START: WRITING TABLE  ===")
		if option.θΨ ≠ "No"
			table.hydroParam.θΨK(Id_Select[1:N_SoilSelect], N_SoilSelect, hydro)
		end
		if option.Psd
			table.psd.PSD(Id_Select[1:N_SoilSelect], N_SoilSelect, psdParam)

			if option.psd.HydroParam
				table.psd.θΨK_PSD(Id_Select, N_SoilSelect, psdHydro)
			end
		end  # if: name
	println("=== END  : WRITING TABLE  === \n")


	if option.Plot
		println("=== START: PLOTTING  ===")

		if option.θΨ ≠ "No" && option.hydro.Plot_θΨ
			plot.HYDROPARAM(Id_Select, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, N_SoilSelect, N_Psd, θ_Rpart, Ψ_Rpart, hydro, psdHydro)
		end # option.Plot_WaterRetentionCurve

		if option.Psd && option.psd.Plot_θr
			plot.PLOT_θr(∑Psd, N_SoilSelect, hydro, psdParam)
		end

		if option.Psd && option.psd.Plot_IMP_model
			plot.PLOT_IMP_model(Id_Select, Rpart, N_Psd, ∑Psd, Psd, N_SoilSelect, hydro, psdParam) 
		end
		# if option.Plot_BestLab && option.θΨ ≠ "No" && (option.infilt.OptimizeRun == "Run" ||  option.infilt.OptimizeRun == "RunOpt") && optiotestn.infiltration.Model=="Simplified" && option.infilt.Plot_SeIni_Range	

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