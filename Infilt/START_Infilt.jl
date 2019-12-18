# =============================================================
#		MODULE: infiltration
# =============================================================
module infilt
	import ..option, ..sorptivity, ..param, ..wrc, ..kunsat, ...opt, ..infiltInitialize, ..bestUniv, ..stats, ..tool
	import BlackBoxOptim, Statistics
	export START_INFILTRATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function START_INFILTRATION(∑Infilt_Obs, ∑Psd, hydro, hydroInfilt, infiltParam, N_Infilt, N_SoilSelect, Tinfilt)

		# INITIALIZE
		T, infiltOutput, hydroInfilt, ∑Infilt = infiltInitialize.INFILT_INITIALIZE(∑Infilt_Obs, ∑Psd, hydroInfilt, infiltParam, N_Infilt, N_SoilSelect, Tinfilt)

		# OPTIONS
			# For every soils
			# hydroInfilt = deepcopy(hydro)
			for iSoil=1:N_SoilSelect
				println( iSoil)

				# No optimization required running from hydro derived from laboratory
				if option.infilt.OptimizeRun == "Run" && option.θΨ ≠ "No" #<>=<>=<>=<>=<>
					# Derive from laboratory
					hydroInfilt = copy(hydro)

					hydroInfilt.θr[iSoil] = min(hydroInfilt.θr[iSoil], infiltParam.θ_Ini[iSoil]) # Not to have errors

					∑Infilt, T_TransStead = INFILTRATION_MODEL(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt, infiltOutput)

				elseif option.infilt.OptimizeRun == "RunOptKs"  && option.θΨ ≠ "No"  #<>=<>=<>=<>=<>
					# Derive from laboratory
				
					SearchRange =[(log10(param.hydro.Ks_Min), log10(param.hydro.Ks_Max))]
					
					hydroInfilt.θr[iSoil] = min(hydroInfilt.θr[iSoil], infiltParam.θ_Ini[iSoil]) # Not to have errors

					Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(∑Infilt, ∑Infilt_Obs, hydroInfilt, infiltParam, iSoil, N_Infilt, T, infiltOutput; Ks=10.0^P[1])[1]; SearchRange=SearchRange, NumDimensions=1, TraceMode=:silent)
		
					hydroInfilt.Ks[iSoil] = 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[1]

					# OUTPUTS
					∑Infilt, T_TransStead = INFILTRATION_MODEL(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt, infiltOutput)

					# Statistics
					 iT_TransSteady = infiltOutput.iT_TransSteady_Data[iSoil]

					infiltOutput.Nse_Trans[iSoil] = 1.0 - stats.NASH_SUTCLIFE_MINIMIZE( ∑Infilt_Obs[iSoil,1:iT_TransSteady], ∑Infilt[iSoil, 1:iT_TransSteady]; Power=2.0)

					infiltOutput.Nse_Steady[iSoil] = 1.0 - stats.NASH_SUTCLIFE_MINIMIZE( log10.(∑Infilt_Obs[iSoil,iT_TransSteady+1:N_Infilt[iSoil]]), log10.(∑Infilt[iSoil,iT_TransSteady+1:N_Infilt[iSoil]]); Power=2.0)

					infiltOutput.Nse[iSoil] = 0.5 * infiltOutput.Nse_Trans[iSoil] + 0.5 * infiltOutput.Nse_Steady[iSoil]

			end # OptimizeRun = "Run"
		end # iSoil

		Nse_Trans = Statistics.mean(infiltOutput.Nse_Trans[1:N_SoilSelect])
		Nse_Steady = Statistics.mean(infiltOutput.Nse_Steady[1:N_SoilSelect])
		Nse =     Statistics.mean(infiltOutput.Nse[1:N_SoilSelect])

		println("    ~ Nse= $Nse Nse_Trans= $Nse_Trans,  Nse_Steady= $Nse_Steady ~")

		println("    ~ Nse= $Nse Nse_Trans= $Nse_Trans,  Nse_Steady= $Nse_Steady ~")

		# CONVERTING DIMENSIONS
			∑Infilt = CONVERT_INFILT_DIMENSIONS(hydroInfilt, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)

		return infiltOutput, hydroInfilt, ∑Infilt
	end  # function: START_INFILTRATION

		
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_MODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_MODEL(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt, infiltOutput)
			if option.infilt.Model == "Best_Univ" # <>=<>=<>=<>=<>

				∑Infilt, T_TransStead = bestUniv.BEST_UNIVERSAL_START(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt, infiltOutput)

			elseif option.infilt.Model == "QuasiExact" # <>=<>=<>=<>=<>
				# quasiExact.QUASIEXACT()

			end #  option.infilt.Model
			
			return ∑Infilt, T_TransStead
		end  # function: INFILTRATION_MODEL


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERT_INFILT_DIMENSIONS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		function CONVERT_INFILT_DIMENSIONS(hydroInfilt, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)
			if option.infilt.Model == "Best_Univ" && option.infilt.SingleDoubleRing == "Double" && option.infilt.OutputDimension == "3D" # <>=<>=<>=<>=<>

				println("    ~ Converting $(option.infilt.Model) Infilt_1D => Infilt_3D ~")

				∑Infilt = bestUniv.CONVERT_1D_2_3D(hydroInfilt, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)

			elseif option.infilt.Model == "Best_Univ" && option.infilt.SingleDoubleRing == "Single" && option.infilt.OutputDimension == "1D" # <>=<>=<>=<>=<>

				println("    ~ Converting $(option.infilt.Model) Infilt_3D => Infilt_1D ~")

				∑Infilt = bestUniv.CONVERT_3D_2_1D(hydroInfilt, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)

			elseif option.infilt.Model == "QuasiExact" # <>=<>=<>=<>=<>
				# quasiExact.QUASIEXACT()
			end #  option.infilt

			return  ∑Infilt
		end  # function: CONVERT_INFILT_DIMENSIONS


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : HYDRO_PROCESS
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OBJECTIVE_FUNCTION(∑Infilt, ∑Infilt_Obs, hydroInfilt, infiltParam, iSoil, N_Infilt, T, infiltOutput; σ=hydroInfilt.σ[iSoil], Ψm=hydroInfilt.Ψm[iSoil], θr=hydroInfilt.θr[iSoil], θs=hydroInfilt.θs[iSoil], Ks=hydroInfilt.Ks[iSoil], W=0.7)

				hydroInfilt.θs[iSoil] = θs
				hydroInfilt.θr[iSoil] = θr
				hydroInfilt.Ks[iSoil] = Ks
				hydroInfilt.σ[iSoil] = σ
				hydroInfilt.Ψm[iSoil] = Ψm

				hydroInfilt.θsMat[iSoil] = θs
				hydroInfilt.ΨmMac[iSoil] = Ψm
				hydroInfilt.σMac[iSoil] = σ

				∑Infilt, T_TransStead = INFILTRATION_MODEL(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt, infiltOutput)

				# Computing the index of T_TransStead
# 				iTime_TransStead = min(sum(T[iSoil,1:N_Infilt[iSoil]] .≤ T_TransStead)+1, N_Infilt[iSoil])

				iT_TransSteady = infiltOutput.iT_TransSteady_Data[iSoil]

				Nse_Trans = stats.NASH_SUTCLIFE_MINIMIZE( ∑Infilt_Obs[iSoil,1:iT_TransSteady], ∑Infilt[iSoil, 1:iT_TransSteady]; Power=2.0)

				Nse_Steady = stats.NASH_SUTCLIFE_MINIMIZE( log10.(∑Infilt_Obs[iSoil,iT_TransSteady+1:N_Infilt[iSoil]]), log10.(∑Infilt[iSoil,iT_TransSteady+1:N_Infilt[iSoil]]); Power=2.0)

				Penalty = abs(∑Infilt_Obs[iSoil,N_Infilt[iSoil]] - ∑Infilt[iSoil,N_Infilt[iSoil] ]) /  ∑Infilt_Obs[iSoil,N_Infilt[iSoil]]

				Nse = W * Nse_Trans + (1.0 - W) * Nse_Steady  + Penalty / 10.0

				return Nse
			end  # function: HYDRO_PROCESS
	

end  # module infiltration
