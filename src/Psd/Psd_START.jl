include("PsdStruct.jl")
include("PsdInitialize.jl")
include("PsdFunc.jl")
include("PsdOpt.jl")


module psdStart
	import ..stats, ..psdInitialize, ..psdOpt, ..psdFunc

	# ======================================================================================
	#          PSD_START Slow initialization
	# ======================================================================================
	function START_PSD(;∑Psd, hydro=[], hydroPsd, NiZ, N_Psd, N_θΨobs, option, param, Rpart, θ_θΨobs, Ψ_θΨobs)

		if  option.psd.OptimizePsd⍰ == "OptAllSoil" && option.run.HydroLabθΨ⍰ == "No"
			error("PSD error: if  option.psd.OptimizePsd⍰ == :OptAllSoil && option.run.HydroLabθΨ⍰ ≠ :No ")
		end 

		# INITIATING THE PSD DATA		
		N_Psd, N_Psd_Max, Psd, θs_Psd, hydroPsd, paramPsd = psdInitialize.PSD_INITIALIZE(∑Psd, hydro, hydroPsd, N_Psd, NiZ, option, param)

		if option.psd.OptimizePsd⍰ == "Run"  # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
			println("\n    == RUN the PSD parameters with prescribed parameters ==") 
			θ_Rpart, Ψ_Rpart = psdOpt.PSD_RUN_ALLMODEL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd,  paramPsd.θr_Psd, option, param, paramPsd, hydroPsd)

			if option.run.HydroLabθΨ⍰ ≠ "No" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				println("\n    == DERIVE HYDAULIC PARAMETERS FROM θ_Rpart, Ψ_Rpart ==")

				paramPsd.Nse, Nse_Mean_Run, Nse_Std_Run = stats.NSE_θΨ(hydroPsd, N_Psd, NiZ, option.psd, θ_Rpart, Ψ_Rpart)
				println("    ~ Nse_Mean_Run=$Nse_Mean_Run, Nse_Std_Run=$Nse_Std_Run ~")
			end
		
			# && option.run.HydroLabθΨ⍰ ≠ :No
		elseif option.psd.OptimizePsd⍰ == "OptSingleSoil" # <>=<>=<>=<>=<>=<>
			if option.psd.Model⍰ == "IMP" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = psdOpt.imp.OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, paramPsd.θr_Psd, option, param, paramPsd, hydro)

			elseif option.psd.Model⍰ == "Chang2019Model" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = psdOpt.chang.OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, paramPsd.θr_Psd, paramPsd, hydro)
			end

			println("\n    == Optimizing the PSD parameters individually for each soil ==")
			println("    	~ Nse_Mean_SingleOpt=$Nse_Mean_SingleOpt,  Nse_Std_SingleOpt=$Nse_Std_SingleOpt ~ \n")

			# && option.run.HydroLabθΨ⍰ ≠ :No
		elseif option.psd.OptimizePsd⍰ == "OptAllSoil" # <>=<>=<>=<>=<>=<>
			if option.psd.Model⍰ == "IMP" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = psdOpt.imp.OPTIMIZATION_ALL_SOIL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, option, θs_Psd, paramPsd.θr_Psd, param, paramPsd, hydro)

			elseif option.psd.Model⍰ == "Chang2019Model" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = psdOpt.chang.OPTIMIZATION_ALL_SOIL(N_Psd_Max, NiZ, Psd, ∑Psd, Rpart, N_Psd, option, θs_Psd, paramPsd.θr_Psd, param, paramPsd, hydro)
			end

			println("\n    == Optimizing the PSD parameters for all soils ==")
			println("    	~ Nse_Mean_OptAllSoil=$Nse_Mean_OptAllSoil,  Nse_Std_OptAllSoil=$Nse_Std_OptAllSoil ~")

  		else
     		error("  $(option.psd.OptimizePsd⍰) not found ")
		end # option.psd.OptimizePsd⍰

	return paramPsd, N_Psd, θ_Rpart, Ψ_Rpart, Psd, hydroPsd	  	
	end # function START_PSD

end # module PSD