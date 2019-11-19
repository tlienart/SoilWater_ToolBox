# =============================================================
#		MODULE: psdOpt
# =============================================================
module psdOpt
	import ..psdFunc

	# =========================================
	#       PSD_RUN_ALLMODEL
	# 		THIS WILL RUN FOR ALL MODELS
	# =========================================
		function PSD_RUN_ALLMODEL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)
			
			θ_Rpart = Array{Float64}(undef, (N_SoilSelect, N_Psd_Max))
			Ψ_Rpart = Array{Float64}(undef, (N_SoilSelect, N_Psd_Max))

			for iSoil = 1:N_SoilSelect
					θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] = psdFunc.PSD_MODEL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam)
			end # for iSoil = 1:N_SoilSelect

			return θ_Rpart, Ψ_Rpart
		end # function PSD_RUN_ALLMODEL

	
	# =============================================================
	#		MODULE: imp
	# =============================================================
	module imp
		import ...param, ...stats, ...wrc, ...option, ...psdFunc, ..psdOpt
		import BlackBoxOptim

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		imp FUNCTION : OPTIMIZATION_SINGLE_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)
				θ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
				Ψ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
				
				for iSoil = 1:N_SoilSelect
					if !(option.psd.∑Psd_2_ξ1)
					
						SearchRange = [ (param.psd.imp.∑Psd_2_ξ2_β1_Min, param.psd.imp.∑Psd_2_ξ2_β1_Max), (param.psd.imp.∑Psd_2_ξ2_β2_Min, param.psd.imp.∑Psd_2_ξ2_β2_Max), (param.psd.imp.Subclay_Min, param.psd.imp.Subclay_Max) ]

						Optimization = BlackBoxOptim.bboptimize(P->OF_SINGLE_SOIL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam, hydro; ∑Psd_2_ξ2_β1 = P[1], ∑Psd_2_ξ2_β2 = P[2], Subclay = P[3])
						; SearchRange = SearchRange, NumDimensions = 3, TraceMode = :silent)

						psdparam.∑Psd_2_ξ2_β1[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
						psdparam.∑Psd_2_ξ2_β2[iSoil] = BlackBoxOptim.best_candidate(Optimization)[2]
						psdparam.Subclay[iSoil]      = BlackBoxOptim.best_candidate(Optimization)[3]


					elseif option.psd.∑Psd_2_ξ1	
						SearchRange =  [(param.psd.imp.ξ1_Min, param.psd.imp.ξ1_Max), (param.psd.imp.∑Psd_2_ξ2_β1_Min, param.psd.imp.∑Psd_2_ξ2_β1_Max), (param.psd.imp.∑Psd_2_ξ2_β2_Min, param.psd.imp.∑Psd_2_ξ2_β2_Max), (param.psd.imp.Subclay_Min, param.psd.imp.Subclay_Max) ]

						Optimization = BlackBoxOptim.bboptimize(P->OF_SINGLE_SOIL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam, hydro; ξ1 = P[1], ∑Psd_2_ξ2_β1 = P[2] ,∑Psd_2_ξ2_β2 = P[3], Subclay = P[4])
						; SearchRange = SearchRange, NumDimensions = 4, TraceMode = :silent)

						psdparam.ξ1[iSoil]           = BlackBoxOptim.best_candidate(Optimization)[1]
						psdparam.∑Psd_2_ξ2_β1[iSoil] = BlackBoxOptim.best_candidate(Optimization)[2]
						psdparam.∑Psd_2_ξ2_β2[iSoil] = BlackBoxOptim.best_candidate(Optimization)[3]
						psdparam.Subclay[iSoil]      = BlackBoxOptim.best_candidate(Optimization)[4]
					end # if option.psd.
				end	# for iSoil = 1:N_SoilSelect
				
				# Compute the optimal values
				θ_Rpart, Ψ_Rpart = psdOpt.PSD_RUN_ALLMODEL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# Statistics
				psdparam.Nse, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydro)

				return psdparam, θ_Rpart, Ψ_Rpart, Nse_Mean_SingleOpt, Nse_Std_SingleOpt
			end # function: OPTIMIZATION_SINGLE_SOIL


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		imp FUNCTION : OPTIMIZATION_ALL_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_ALL_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				#		imp FUNCTION : OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					function OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, N_SoilSelect, θs_Psd, θr_Psd, psdparam, hydro; ξ1 = param.psd.imp.ξ1, ∑Psd_2_ξ2_β1 = param.psd.imp.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2 = param.psd.imp.∑Psd_2_ξ2_β2, Subclay = param.psd.imp.Subclay)

						θ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
						Ψ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))

						for iSoil = 1:N_SoilSelect
							psdparam.ξ1[iSoil]           = ξ1
							psdparam.∑Psd_2_ξ2_β1[iSoil] = ∑Psd_2_ξ2_β1
							psdparam.∑Psd_2_ξ2_β2[iSoil] = ∑Psd_2_ξ2_β2
							psdparam.Subclay[iSoil]      = Subclay
						end

						Of = 0.0
						for iSoil = 1:N_SoilSelect
							θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] = psdFunc.PSD_MODEL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam)

							θΨ =  Array{Float64}(undef, (N_Psd[iSoil]))
							for iRpart = 1:N_Psd[iSoil]
								# Observed data
								θΨ[iRpart] = wrc.Ψ_2_θDual(Ψ_Rpart[iSoil, iRpart], iSoil, hydro)
							end

							Of += stats.NASH_SUTCLIFE_MINIMIZE(θΨ[1:N_Psd[iSoil]], θ_Rpart[iSoil,1:N_Psd[iSoil]]; Power=2)
						end # for iSoil = 1:N_SoilSelect
						return Of
					end  # function OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


				if !(option.psd.∑Psd_2_ξ1)	
					SearchRange = [ (param.psd.imp.∑Psd_2_ξ2_β1_Min, param.psd.imp.∑Psd_2_ξ2_β1_Max), (param.psd.imp.∑Psd_2_ξ2_β2_Min, param.psd.imp.∑Psd_2_ξ2_β2_Max), (param.psd.imp.Subclay_Min, param.psd.imp.Subclay_Max) ]

					Optimization = BlackBoxOptim.bboptimize(P->OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, N_SoilSelect, θs_Psd, θr_Psd, psdparam, hydro; ∑Psd_2_ξ2_β1 = P[1], ∑Psd_2_ξ2_β2 = P[2], Subclay = P[3])
					; SearchRange = SearchRange, NumDimensions = 3, TraceMode = :silent)

					for iSoil=N_SoilSelect
						psdparam.∑Psd_2_ξ2_β1[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
						psdparam.∑Psd_2_ξ2_β2[iSoil] = BlackBoxOptim.best_candidate(Optimization)[2]
						psdparam.Subclay[iSoil]      = BlackBoxOptim.best_candidate(Optimization)[3]
					end

				elseif option.psd.∑Psd_2_ξ1
					SearchRange =  [(param.psd.imp.ξ1_Min, param.psd.imp.ξ1_Max), (param.psd.imp.∑Psd_2_ξ2_β1_Min, param.psd.imp.∑Psd_2_ξ2_β1_Max), (param.psd.imp.∑Psd_2_ξ2_β2_Min, param.psd.imp.∑Psd_2_ξ2_β2_Max), (param.psd.imp.Subclay_Min, param.psd.imp.Subclay_Max) ]

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
				θ_Rpart, Ψ_Rpart = psdOpt.PSD_RUN_ALLMODEL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# Statistics
				psdparam.Nse, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydro)

				return psdparam, θ_Rpart, Ψ_Rpart, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil
			end # function OPTIMIZATION_ALL_SOIL


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		imp FUNCTION : OF_SINGLE_SOIL
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function OF_SINGLE_SOIL(iSoil, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro; ξ1 = psdparam.ξ1[iSoil], ∑Psd_2_ξ2_β1 = psdparam.∑Psd_2_ξ2_β1[iSoil], ∑Psd_2_ξ2_β2 = psdparam.∑Psd_2_ξ2_β2[iSoil], Subclay = psdparam.Subclay[iSoil])

					psdparam.ξ1[iSoil]           = ξ1
					psdparam.∑Psd_2_ξ2_β1[iSoil] = ∑Psd_2_ξ2_β1
					psdparam.∑Psd_2_ξ2_β2[iSoil] = ∑Psd_2_ξ2_β2
					psdparam.Subclay[iSoil]      = Subclay

					# Compute the proposed value
					θ_Rpart, Ψ_Rpart = psdFunc.PSD_MODEL(iSoil, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam)

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
	 	import ...param, ...stats, ...wrc, ...option, ....psdFunc, ..psdOpt
		import BlackBoxOptim

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		chang FUNCTION : OPTIMIZATION_SINGLE_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)
				θ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
				Ψ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
				
				for iSoil = 1:N_SoilSelect
					SearchRange = [ (param.psd.chan.ξ1_Min, param.psd.chan.ξ1_Max) ]

					Optimization = BlackBoxOptim.bboptimize(P->OF_SINGLE_SOIL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam, hydro; ξ1=P[1])
					; SearchRange = SearchRange, NumDimensions=1, TraceMode = :silent)

					psdparam.ξ1[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
				end	# for iSoil = 1:N_SoilSelect
				
				# **Compute the optimal values**
				θ_Rpart, Ψ_Rpart = psdOpt.PSD_RUN_ALLMODEL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# Statistics
				psdparam.Nse, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydro)

				return psdparam, θ_Rpart, Ψ_Rpart, Nse_Mean_SingleOpt, Nse_Std_SingleOpt
			end # function: OPTIMIZATION_SINGLE_SOIL


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		chang FUNCTION : OPTIMIZATION_ALL_SOIL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMIZATION_ALL_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				#		chang FUNCTION : OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					function OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, N_SoilSelect, θs_Psd, θr_Psd, psdparam, hydro; ξ1 = psdparam.ξ1)
						θ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))
						Ψ_Rpart = zeros(Float64, (N_SoilSelect, N_Psd_Max))

						for iSoil = 1:N_SoilSelect
							psdparam.ξ1[iSoil] = ξ1
						end

						Of = 0.0
						for iSoil = 1:N_SoilSelect
							θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] = psdFunc.PSD_MODEL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θs_Psd[iSoil], θr_Psd[iSoil], psdparam)

							θΨ = zeros(Float64, N_Psd[iSoil])
							for iRpart = 1:N_Psd[iSoil]
								# Observed data
								θΨ[iRpart] = wrc.Ψ_2_θDual(Ψ_Rpart[iSoil, iRpart], iSoil, hydro)
							end

							Of += stats.NASH_SUTCLIFE_MINIMIZE(θΨ[1:N_Psd[iSoil]], θ_Rpart[iSoil,1:N_Psd[iSoil]]; Power=2)
						end # for iSoil = 1:N_SoilSelect
						return Of
					end  # function OF_ALL_SOIL
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				SearchRange = [ (param.psd.chan.ξ1_Min, param.psd.chan.ξ1_Max) ]

				Optimization = BlackBoxOptim.bboptimize(P->OF_ALL_SOIL(Psd, ∑Psd, Rpart, N_Psd, N_Psd_Max, N_SoilSelect, θs_Psd, θr_Psd, psdparam, hydro; ξ1 = P[1])
				; SearchRange = SearchRange, NumDimensions=1, TraceMode = :silent)

				for iSoil=1:N_SoilSelect
					psdparam.ξ1[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
				end

				# **Compute the optimal values**
				θ_Rpart, Ψ_Rpart = psdOpt.PSD_RUN_ALLMODEL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro)

				# Statistics
				psdparam.Nse, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydro)

				return psdparam, θ_Rpart, Ψ_Rpart, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil
			end # function OPTIMIZATION_ALL_SOIL


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		chang FUNCTION : OF_SINGLE_SOIL
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function OF_SINGLE_SOIL(iSoil, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam, hydro; ξ1 = psdparam.ξ1[iSoil])

					psdparam.ξ1[iSoil] = ξ1

					# Compute the proposed value
					θ_Rpart, Ψ_Rpart = psdFunc.PSD_MODEL(iSoil, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, psdparam)

					# For every Ψ_Rpart
					θΨ = zeros(Float64, N_Psd)
					for iRpart = 1:N_Psd
						# Observed data
						θΨ[iRpart] = wrc.Ψ_2_θDual(Ψ_Rpart[iRpart], iSoil, hydro)
					end

					Of = stats.NASH_SUTCLIFE_MINIMIZE(θΨ[1:N_Psd], θ_Rpart[1:N_Psd]; Power=2)
				end  # function OF_SINGLE_SOIL
		
		
	end  # module chang
	# ............................................................


end  # module psdOpt
# ............................................................