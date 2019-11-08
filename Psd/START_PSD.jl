module psd
	import ..option, ..psdInitiate, ..psdThetar, ..psdStruct

	# ======================================================================================
	#          PSD_START
	# ======================================================================================
	function START_PSD(N_SoilSelect, Ψ_θΨ, θ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, Rpart, ∑Psd, N_Psd, hydro)
		# INITIATING THE PSD DATA		
 			N_Psd, N_Psd_Max, Psd = psdInitiate.PSD_INITIATE(N_Psd, N_SoilSelect, ∑Psd)
		
		# COMPUTING θr FROM PSD DATA
 			Nse_θr, θr_Psd = psdThetar.MAIN_PSD_2_θr(N_SoilSelect, ∑Psd, hydro)

 			θs_Psd = hydro.θs[1:N_SoilSelect] # TODO need to read θs from file

		# DERIVING THE STRUCTURE PARAMETERS
			 psdparam = psdStruct.PSDSTRUCT(N_SoilSelect)

		if option.psd.OptimizePsd == "Run"  # <> <> <> <> <> <> 
			θ_Rpart, Ψ_Rpart = allmodel.PSD_RUN(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)
			
			if option.θΨ ≠ "No"
				Nse_Run, Nse_Mean_Run, Nse_Std_Run = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydro)
				println("    ~ Nse_Mean_Run=$Nse_Mean_Run, Nse_Std_Run=$Nse_Std_Run ~")
			end
		

		elseif option.psd.OptimizePsd == "OptSingleSoil" && option.θΨ ≠ "No" # <> <> <> <> <> <>
			psdparam, θ_Rpart, Ψ_Rpart, Nse_SingleOpt, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = imp.OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

			println("    ~ Nse_Mean_SingleOpt=$Nse_Mean_SingleOpt,  Nse_Std_SingleOpt=$Nse_Std_SingleOpt ~")


		elseif option.psd.OptimizePsd == "OptAllSoil" && option.θΨ ≠ "No" # <> <> <> <> <> <>
			psdparam, θ_Rpart, Ψ_Rpart, Nse_OptAllSoil, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = imp.OPTIMIZATION_ALL_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

			println("    ~ Nse_Mean_OptAllSoil=$Nse_Mean_OptAllSoil,  Nse_Std_OptAllSoil=$Nse_Std_OptAllSoil ~")

  		else
     		error("  $(option.psd.OptimizePsd) not found ")
		end # option.psd.OptimizePsd

		return psdparam, θ_Rpart, Ψ_Rpart
		  	
	end # function START_PSD

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# =============================================================
	#		MODULE: allmodels
	# =============================================================
	module allmodel
		import  ...psdFunc

		# =========================================
		#       PSD_RUN
		# 		THIS WILL RUN FOR ALL MODELS
		# =========================================
		function PSD_RUN(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)
			θ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
			Ψ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))

			for iSoil = 1:N_SoilSelect
				   θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] = psdFunc._PSD_MODEL_(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam)
			end # for iSoil = 1:N_SoilSelect

			return θ_Rpart, Ψ_Rpart
		end # function PSD_RUN

	end  # module allmodels
	# ............................................................

	# =============================================================
	#		MODULE: imp
	# =============================================================
	module imp
		import ...param, ...stats, ...wrc, ...psdFunc, ..option, ..allmodel
		import BlackBoxOptim

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : OPTIMIZATION_SINGLE_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)
				θ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
				Ψ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
				
				@simd for iSoil = 1:N_SoilSelect
					if !(option.psd.∑Psd_2_ξ1)
					
						SearchRange = [ (param.psd.∑Psd_2_ξ2_β1_Min, param.psd.∑Psd_2_ξ2_β1_Max), (param.psd.∑Psd_2_ξ2_β2_Min, param.psd.∑Psd_2_ξ2_β2_Max), (param.psd.Subclay_Min, param.psd.Subclay_Max) ]

						Optimization = BlackBoxOptim.bboptimize(P->OF_SINGLE_SOIL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam, hydro; ∑Psd_2_ξ2_β1 = P[1], ∑Psd_2_ξ2_β2 = P[2], Subclay = P[3])
						; SearchRange = SearchRange, NumDimensions = 3, TraceMode = :silent)

						psdparam.∑Psd_2_ξ2_β1[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
						psdparam.∑Psd_2_ξ2_β2[iSoil] = BlackBoxOptim.best_candidate(Optimization)[2]
						psdparam.Subclay[iSoil]      = BlackBoxOptim.best_candidate(Optimization)[3]


					elseif option.psd.∑Psd_2_ξ1	
						SearchRange =  [(param.psd.ξ1_Min, param.psd.ξ1_Max), (param.psd.∑Psd_2_ξ2_β1_Min, param.psd.∑Psd_2_ξ2_β1_Max), (param.psd.∑Psd_2_ξ2_β2_Min, param.psd.∑Psd_2_ξ2_β2_Max), (param.psd.Subclay_Min, param.psd.Subclay_Max) ]

						Optimization = BlackBoxOptim.bboptimize(P->OF_SINGLE_SOIL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam, hydro; ξ1 = P[1], ∑Psd_2_ξ2_β1 = P[2] ,∑Psd_2_ξ2_β2 = P[3], Subclay = P[4])
						; SearchRange = SearchRange, NumDimensions = 4, TraceMode = :silent)

						psdparam.ξ1[iSoil]           = BlackBoxOptim.best_candidate(Optimization)[1]
						psdparam.∑Psd_2_ξ2_β1[iSoil] = BlackBoxOptim.best_candidate(Optimization)[2]
						psdparam.∑Psd_2_ξ2_β2[iSoil] = BlackBoxOptim.best_candidate(Optimization)[3]
						psdparam.Subclay[iSoil]      = BlackBoxOptim.best_candidate(Optimization)[4]
					end # if option.psd.
				end	# for iSoil = 1:N_SoilSelect
				
				# **Compute the optimal values**
				θ_Rpart, Ψ_Rpart = PSD_RUN(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# Statistics
				Nse_SingleOpt, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydro)

				return psdparam, θ_Rpart, Ψ_Rpart, Nse_SingleOpt, Nse_Mean_SingleOpt, Nse_Std_SingleOpt
			end # function: OPTIMIZATION_SINGLE_SOIL


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : OPTIMIZATION_ALL_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_ALL_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				#		FUNCTION : OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					function OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, N_SoilSelect, θs_Psd, θr_Psd, psdparam, hydro; ξ1 = psdparam.ξ1, ∑Psd_2_ξ2_β1 = psdparam.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2 = psdparam.∑Psd_2_ξ2_β2, Subclay = psdparam.Subclay)
						θ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
						Ψ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))

						@simd for iSoil = 1:N_SoilSelect
							psdparam.ξ1[iSoil]           = ξ1
							psdparam.∑Psd_2_ξ2_β1[iSoil] = ∑Psd_2_ξ2_β1
							psdparam.∑Psd_2_ξ2_β2[iSoil] = ∑Psd_2_ξ2_β2
							psdparam.Subclay[iSoil]      = Subclay
						end

						Of = 0.0
						@simd for iSoil = 1:N_SoilSelect
							θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] = psdFunc._PSD_MODEL_(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam)

							θΨ = zeros(Float64, N_Psd[iSoil])
							for iRpart = 1:N_Psd[iSoil]
								# Observed data
								θΨ[iRpart] = wrc.Ψ_2_θDual(Ψ_Rpart[iSoil, iRpart], iSoil, hydro)
							end

							Of += stats.NASH_SUTCLIFE_MINIMIZE(θΨ[1:N_Psd[iSoil]], θ_Rpart[1:N_Psd[iSoil]]; Power=2)
						end # for iSoil = 1:N_SoilSelect
						return Of
					end  # function OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


				if !(option.psd.∑Psd_2_ξ1)	
					SearchRange = [ (param.psd.∑Psd_2_ξ2_β1_Min, param.psd.∑Psd_2_ξ2_β1_Max), (param.psd.∑Psd_2_ξ2_β2_Min, param.psd.∑Psd_2_ξ2_β2_Max), (param.psd.Subclay_Min, param.psd.Subclay_Max) ]

					Optimization = BlackBoxOptim.bboptimize(P->OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, N_SoilSelect, θs_Psd, θr_Psd, psdparam, hydro; ∑Psd_2_ξ2_β1 = P[1], ∑Psd_2_ξ2_β2 = P[2], Subclay = P[3])
					; SearchRange = SearchRange, NumDimensions = 3, TraceMode = :silent)

					for iSoil=N_SoilSelect
						psdparam.∑Psd_2_ξ2_β1[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
						psdparam.∑Psd_2_ξ2_β2[iSoil] = BlackBoxOptim.best_candidate(Optimization)[2]
						psdparam.Subclay[iSoil]      = BlackBoxOptim.best_candidate(Optimization)[3]
					end

				elseif option.psd.∑Psd_2_ξ1
					SearchRange =  [(param.psd.ξ1_Min, param.psd.ξ1_Max), (param.psd.∑Psd_2_ξ2_β1_Min, param.psd.∑Psd_2_ξ2_β1_Max), (param.psd.∑Psd_2_ξ2_β2_Min, param.psd.∑Psd_2_ξ2_β2_Max), (param.psd.Subclay_Min, param.psd.Subclay_Max) ]

					Optimization = BlackBoxOptim.bboptimize(P->OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, N_SoilSelect, θs_Psd, θr_Psd, psdparam, hydro; ξ1 = P[1], ∑Psd_2_ξ2_β1 = P[2] ,∑Psd_2_ξ2_β2 = P[3], Subclay = P[4])
					; SearchRange = SearchRange, NumDimensions = 4, TraceMode = :silent)

					for iSoil=N_SoilSelect
						psdparam.ξ1[iSoil]           = BlackBoxOptim.best_candidate(Optimization)[1]
						psdparam.∑Psd_2_ξ2_β1[iSoil] = BlackBoxOptim.best_candidate(Optimization)[2]
						psdparam.∑Psd_2_ξ2_β2[iSoil] = BlackBoxOptim.best_candidate(Optimization)[3]
						psdparam.Subclay[iSoil]      = BlackBoxOptim.best_candidate(Optimization)[4]
					end
				end # if option.psd.

				# **Compute the optimal values**
				θ_Rpart, Ψ_Rpart = allmodel.PSD_RUN(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# Statistics
				Nse_OptAllSoil, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydro)

				return psdparam, θ_Rpart, Ψ_Rpart, Nse_OptAllSoil, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil
			end # function OPTIMIZATION_ALL_SOIL


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : OF_SINGLE_SOIL
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function OF_SINGLE_SOIL(iSoil, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro; ξ1 = psdparam.ξ1[iSoil], ∑Psd_2_ξ2_β1 = psdparam.∑Psd_2_ξ2_β1[iSoil], ∑Psd_2_ξ2_β2 = psdparam.∑Psd_2_ξ2_β2[iSoil], Subclay = psdparam.Subclay[iSoil])

					psdparam.ξ1[iSoil]           = ξ1
					psdparam.∑Psd_2_ξ2_β1[iSoil] = ∑Psd_2_ξ2_β1
					psdparam.∑Psd_2_ξ2_β2[iSoil] = ∑Psd_2_ξ2_β2
					psdparam.Subclay[iSoil]      = Subclay

					# Compute the proposed value
					θ_Rpart, Ψ_Rpart = psdFunc._PSD_MODEL_(iSoil, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam)

					# For every Ψ_Rpart
					Of = 0.0
					θΨ = zeros(Float64, N_Psd)
					for iRpart = 1:N_Psd
						# Observed data
						θΨ[iRpart] = wrc.Ψ_2_θDual(Ψ_Rpart[iRpart], iSoil, hydro)
					end

					Of = stats.NASH_SUTCLIFE_MINIMIZE(θΨ[1:N_Psd], θ_Rpart[1:N_Psd]; Power=2)
				end  # function OF_SINGLE_SOIL
		
	end  # module imp
	# ............................................................


	# =============================================================
	#		MODULE: chang
	# =============================================================
	module chang

	import ...param, ...stats, ...wrc, ...psdFunc, ..option, ..allmodel
		import BlackBoxOptim

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : OPTIMIZATION_SINGLE_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)
				θ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
				Ψ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
				
				@simd for iSoil = 1:N_SoilSelect
					SearchRange = [ (param.psd.∑Psd_2_ξ2_β1_Min, param.psd.∑Psd_2_ξ2_β1_Max), (param.psd.∑Psd_2_ξ2_β2_Min, param.psd.∑Psd_2_ξ2_β2_Max), (param.psd.Subclay_Min, param.psd.Subclay_Max) ]

					Optimization = BlackBoxOptim.bboptimize(P->OF_SINGLE_SOIL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam, hydro; ∑Psd_2_ξ2_β1 = P[1], ∑Psd_2_ξ2_β2 = P[2], Subclay = P[3])
					; SearchRange = SearchRange, NumDimensions = 3, TraceMode = :silent)

					psdparam.∑Psd_2_ξ2_β1[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
					psdparam.∑Psd_2_ξ2_β2[iSoil] = BlackBoxOptim.best_candidate(Optimization)[2]
					psdparam.Subclay[iSoil]      = BlackBoxOptim.best_candidate(Optimization)[3]

				end	# for iSoil = 1:N_SoilSelect
				
				# **Compute the optimal values**
				θ_Rpart, Ψ_Rpart = PSD_RUN(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# Statistics
				Nse_SingleOpt, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydro)

				return psdparam, θ_Rpart, Ψ_Rpart, Nse_SingleOpt, Nse_Mean_SingleOpt, Nse_Std_SingleOpt
			end # function: OPTIMIZATION_SINGLE_SOIL


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : OPTIMIZATION_ALL_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_ALL_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				#		FUNCTION : OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					function OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, N_SoilSelect, θs_Psd, θr_Psd, psdparam, hydro; ξ1 = psdparam.ξ1, ∑Psd_2_ξ2_β1 = psdparam.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2 = psdparam.∑Psd_2_ξ2_β2, Subclay = psdparam.Subclay)
						θ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
						Ψ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))

						@simd for iSoil = 1:N_SoilSelect
							psdparam.ξ1[iSoil]           = ξ1
							psdparam.∑Psd_2_ξ2_β1[iSoil] = ∑Psd_2_ξ2_β1
							psdparam.∑Psd_2_ξ2_β2[iSoil] = ∑Psd_2_ξ2_β2
							psdparam.Subclay[iSoil]      = Subclay
						end

						Of = 0.0
						@simd for iSoil = 1:N_SoilSelect
							θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] = psdFunc._PSD_MODEL_(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam)

							θΨ = zeros(Float64, N_Psd[iSoil])
							for iRpart = 1:N_Psd[iSoil]
								# Observed data
								θΨ[iRpart] = wrc.Ψ_2_θDual(Ψ_Rpart[iSoil, iRpart], iSoil, hydro)
							end

							Of += stats.NASH_SUTCLIFE_MINIMIZE(θΨ[1:N_Psd[iSoil]], θ_Rpart[1:N_Psd[iSoil]]; Power=2)
						end # for iSoil = 1:N_SoilSelect
						return Of
					end  # function OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				SearchRange = [ (param.psd.∑Psd_2_ξ2_β1_Min, param.psd.∑Psd_2_ξ2_β1_Max), (param.psd.∑Psd_2_ξ2_β2_Min, param.psd.∑Psd_2_ξ2_β2_Max), (param.psd.Subclay_Min, param.psd.Subclay_Max) ]

				Optimization = BlackBoxOptim.bboptimize(P->OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, N_SoilSelect, θs_Psd, θr_Psd, psdparam, hydro; ∑Psd_2_ξ2_β1 = P[1], ∑Psd_2_ξ2_β2 = P[2], Subclay = P[3])
				; SearchRange = SearchRange, NumDimensions = 3, TraceMode = :silent)

				for iSoil=N_SoilSelect
					psdparam.∑Psd_2_ξ2_β1[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
					psdparam.∑Psd_2_ξ2_β2[iSoil] = BlackBoxOptim.best_candidate(Optimization)[2]
					psdparam.Subclay[iSoil]      = BlackBoxOptim.best_candidate(Optimization)[3]
				end

			

				# **Compute the optimal values**
				θ_Rpart, Ψ_Rpart = allmodel.PSD_RUN(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# Statistics
				Nse_OptAllSoil, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydro)

				return psdparam, θ_Rpart, Ψ_Rpart, Nse_OptAllSoil, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil
			end # function OPTIMIZATION_ALL_SOIL


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : OF_SINGLE_SOIL
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function OF_SINGLE_SOIL(iSoil, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro; ξ1 = psdparam.ξ1[iSoil], ∑Psd_2_ξ2_β1 = psdparam.∑Psd_2_ξ2_β1[iSoil], ∑Psd_2_ξ2_β2 = psdparam.∑Psd_2_ξ2_β2[iSoil], Subclay = psdparam.Subclay[iSoil])

					psdparam.ξ1[iSoil]           = ξ1
					psdparam.∑Psd_2_ξ2_β1[iSoil] = ∑Psd_2_ξ2_β1
					psdparam.∑Psd_2_ξ2_β2[iSoil] = ∑Psd_2_ξ2_β2
					psdparam.Subclay[iSoil]      = Subclay

					# Compute the proposed value
					θ_Rpart, Ψ_Rpart = psdFunc._PSD_MODEL_(iSoil, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam)

					# For every Ψ_Rpart
					Of = 0.0
					θΨ = zeros(Float64, N_Psd)
					for iRpart = 1:N_Psd
						# Observed data
						θΨ[iRpart] = wrc.Ψ_2_θDual(Ψ_Rpart[iRpart], iSoil, hydro)
					end

					Of = stats.NASH_SUTCLIFE_MINIMIZE(θΨ[1:N_Psd], θ_Rpart[1:N_Psd]; Power=2)
				end  # function OF_SINGLE_SOIL
		
		
	end  # module chang
	# ............................................................

end # module PSD