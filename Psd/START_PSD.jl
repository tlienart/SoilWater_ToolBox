module psd
	import ..option, ..psdInitialize, ..psdThetar, ..psdStruct, ..stats, ..psdOpt, ..psdFunc, ..param

	# ======================================================================================
	#          PSD_START Slow initialization
	# ======================================================================================
	function START_PSD(N_SoilSelect, Ψ_θΨ, θ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, Rpart, ∑Psd, N_Psd, Φ_Psd, hydro)
		# INITIATING THE PSD DATA		
 			N_Psd, N_Psd_Max, Psd = psdInitialize.PSD_INITIALIZE(N_Psd, N_SoilSelect, ∑Psd)
		
		# COMPUTING θr FROM PSD DATA
 			Nse_θr, θr_Psd = psdThetar.MAIN_PSD_2_θr(N_SoilSelect, ∑Psd, hydro)

 			θs_Psd = param.hydro.Coeff_Φ_2_θs .* Φ_Psd[1:N_SoilSelect]  

		# DERIVING THE STRUCTURE PARAMETERS
			 psdparam = psdStruct.PSDSTRUCT(N_SoilSelect)

		if option.psd.OptimizePsd == "Run"  # <> <> <> <> <> <> 
			θ_Rpart, Ψ_Rpart = psd.PSD_RUN_ALLMODEL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)
			
			if option.θΨ ≠ "No"
				Nse_Run, Nse_Mean_Run, Nse_Std_Run = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydro)
				println("    ~ Nse_Mean_Run=$Nse_Mean_Run, Nse_Std_Run=$Nse_Std_Run ~")
			end
		

		elseif option.psd.OptimizePsd == "OptSingleSoil" && option.θΨ ≠ "No" # <> <> <> <> <> <>
			if option.psd.Model == "IMP"
				psdparam, θ_Rpart, Ψ_Rpart, Nse_SingleOpt, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = psdOpt.imp.OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

			elseif option.psd.Model == "Chang2019Model"
				psdparam, θ_Rpart, Ψ_Rpart, Nse_SingleOpt, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = psdOpt.chang.OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)
			end

			println("    ~ Nse_Mean_SingleOpt=$Nse_Mean_SingleOpt,  Nse_Std_SingleOpt=$Nse_Std_SingleOpt ~")


		elseif option.psd.OptimizePsd == "OptAllSoil" && option.θΨ ≠ "No" # <> <> <> <> <> <>
			if option.psd.Model == "IMP"
				psdparam, θ_Rpart, Ψ_Rpart, Nse_OptAllSoil, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = psdOpt.imp.OPTIMIZATION_ALL_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

			elseif option.psd.Model == "Chang2019Model"
				psdparam, θ_Rpart, Ψ_Rpart, Nse_OptAllSoil, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = psdOpt.chang.OPTIMIZATION_ALL_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)
			end

			println("    ~ Nse_Mean_OptAllSoil=$Nse_Mean_OptAllSoil,  Nse_Std_OptAllSoil=$Nse_Std_OptAllSoil ~")

  		else
     		error("  $(option.psd.OptimizePsd) not found ")
		end # option.psd.OptimizePsd

		return psdparam, θ_Rpart, Ψ_Rpart
		  	
	end # function START_PSD

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>


	# =========================================
	#       PSD_RUN_ALLMODEL
	# 		THIS WILL RUN FOR ALL MODELS
	# =========================================
		function PSD_RUN_ALLMODEL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)
			θ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
			Ψ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))

			for iSoil = 1:N_SoilSelect
					θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] = psdFunc.PSD_MODEL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam)
			end # for iSoil = 1:N_SoilSelect

			return θ_Rpart, Ψ_Rpart
		end # function PSD_RUN_ALLMODEL

end # module PSD