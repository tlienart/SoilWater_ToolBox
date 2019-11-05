module psd
	import ..option, ..param, ..wrc, ..kunsat, ..cst, ..path, ..stats, ..psdFunc, ..psdInitiate, ..psdThetar, ..psdStruct
	using Statistics, BlackBoxOptim, JuliaDB

	# ======================================================================================
	#          PSD_START
	# ======================================================================================
	function START_PSD(N_SoilSelect, Ψ_θΨ, θ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, Rpart, ∑Psd, N_Psd, hydro)

		# INITIATING THE PSD DATA		
			N_Psd, N_Psd_Max, Psd = psdInitiate.PSD_INITIATE(N_Psd, N_SoilSelect, ∑Psd)
		
		# COMPUTING θr FROM PSD DATA
			Nse_θr, θr_Psd = psdThetar.MAIN_PSD_2_θr(N_SoilSelect, ∑Psd, hydro)
			println("NASH_SUTCLIFFE θr_Psd = $Nse_θr \n")

			θs_Psd = hydro.θs[1:N_SoilSelect] # TODO need to read θs from file

		# DERRIVING THE STRUCTURE PARAMETERS
			psdparam = psdStruct.PSDSTRUCT(N_SoilSelect)


		if option.psd.OptimizePsd == "Run" # 

			θ_Rpart, Ψ_Rpart = PSD_RUN(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam)
		
		# if option.psd.OptimizePsd == "Single" # <> <> <> <> <> <>
		# 	∑Of_Psd = 0.0

		# 	@simd for iSoil=1:N_SoilSelect
		# 		ξ1[iSoil], ξ2[iSoil], Ψ_Rpart[iSoil,1:N_Psd[iSoil]], θ_Rpart[iSoil,1:N_Psd[iSoil]], Of_Psd[iSoil], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Subclay[iSoil] = OPTIMIZATION_SINGLE_SOIL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], hydro, θr_Psd[iSoil])
		# 	end # Loop single
``
		# elseif option.psd.OptimizePsd == "All" # <> <> <> <> <> <>
		# 	∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd = OPTIMIZATION_ALL_SOIL(N_SoilSelect, Psd[1:N_SoilSelect,:], ∑Psd[1:N_SoilSelect,:], Rpart[1:N_SoilSelect,:], N_Psd[1:N_SoilSelect], hydro, θr_Psd[1:N_SoilSelect])

		else
			error("  $(option.psd.OptimizePsd) not found ")
		end # option.psd.OptimizePsd	
	end # function START_PSD

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

		# =========================================
		#       PSD_RUN
		# 		THIS WILL RUN FOR ALL MODELS
		# =========================================
			function PSD_RUN( N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam)
				θ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
				Ψ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))

				@simd for iSoil = 1:N_SoilSelect
					θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] =  psdFunc.PSD_MODEL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam)

					println(θ_Rpart[iSoil,1:N_Psd[iSoil]])
					println( Ψ_Rpart[iSoil,1:N_Psd[iSoil]])

				end # for iSoil = 1:N_SoilSelect

				return θ_Rpart, Ψ_Rpart
			end # function PSD_RUN




end # module PSD