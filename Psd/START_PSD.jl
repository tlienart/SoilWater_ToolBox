module psd
	import ..option, ..param, ..stats, ..psdInitialize, ..psdOpt, ..psdFunc

	# ======================================================================================
	#          PSD_START Slow initialization
	# ======================================================================================
	function START_PSD(N_SoilSelect, Ψ_θΨ, θ_θΨ, N_θΨ, Rpart, ∑Psd, N_Psd, Φ_Psd, hydro, psdHydro)

		# INITIATING THE PSD DATA		
		N_Psd, N_Psd_Max, Psd, θs_Psd, psdHydro, psdParam = psdInitialize.PSD_INITIALIZE(N_Psd, N_SoilSelect, ∑Psd, Φ_Psd, hydro, psdHydro)

		if option.psd.OptimizePsd == "Run"  # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
			θ_Rpart, Ψ_Rpart = psdOpt.PSD_RUN_ALLMODEL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, psdParam.θr_Psd, psdParam, hydro)
			
			if option.θΨ ≠ "No" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				psdParam.Nse, Nse_Mean_Run, Nse_Std_Run = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydro)
				println("    ~ Nse_Mean_Run=$Nse_Mean_Run, Nse_Std_Run=$Nse_Std_Run ~")
			end
		

		elseif option.psd.OptimizePsd == "OptSingleSoil" && option.θΨ ≠ "No" # <> <> <> <> <> <>
			if option.psd.Model == "IMP" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				psdParam, θ_Rpart, Ψ_Rpart, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = psdOpt.imp.OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, psdParam.θr_Psd, psdParam, hydro)

			elseif option.psd.Model == "Chang2019Model" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				psdParam, θ_Rpart, Ψ_Rpart, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = psdOpt.chang.OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, psdParam.θr_Psd, psdParam, hydro)
			end

			println("    ~ Nse_Mean_SingleOpt=$Nse_Mean_SingleOpt,  Nse_Std_SingleOpt=$Nse_Std_SingleOpt ~")


		elseif option.psd.OptimizePsd == "OptAllSoil" && option.θΨ ≠ "No" # <> <> <> <> <> <>
			if option.psd.Model == "IMP" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				psdParam, θ_Rpart, Ψ_Rpart, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = psdOpt.imp.OPTIMIZATION_ALL_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, psdParam.θr_Psd, psdParam, hydro)

			elseif option.psd.Model == "Chang2019Model" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				psdParam, θ_Rpart, Ψ_Rpart, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = psdOpt.chang.OPTIMIZATION_ALL_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, psdParam.θr_Psd, psdParam, hydro)
			end

			println("    ~ Nse_Mean_OptAllSoil=$Nse_Mean_OptAllSoil,  Nse_Std_OptAllSoil=$Nse_Std_OptAllSoil ~")

  		else
     		error("  $(option.psd.OptimizePsd) not found ")
		end # option.psd.OptimizePsd

		return psdParam, N_Psd, θ_Rpart, Ψ_Rpart, Psd, psdHydro
		  	
	end # function START_PSD


end # module PSD