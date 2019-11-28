# =============================================================
#		MODULE: infiltration
# =============================================================
module infilt
	import ..option, ..sorptivity, ..param, ..wrc, ..kunsat, ...opt, ..infiltInitialize, ..bestUniv
	import BlackBoxOptim
	export START_INFILTRATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function START_INFILTRATION(N_SoilSelect, Tinfilt, ∑Infilt_Obs, ∑Psd, N_Infilt, infiltParam, hydro)

		# INITIALIZE
			∑Infilt, T, hydroInfilt = infiltInitialize.INFILT_INITIALIZE(N_SoilSelect, ∑Psd, infiltParam, Tinfilt, N_Infilt)

		# OPTIONS
			# For every soils
			hydroInfilt = hydro
			for iSoil=1:N_SoilSelect

				# No optimization required running from hydro derived from laboratory
				if option.infilt.OptimizeRun == "Run" && option.θΨ ≠ "No" #<>=<>=<>=<>=<>
					# Derive from laboratory
					hydroInfilt = hydro

					hydroInfilt.θr[iSoil] = min(hydroInfilt.θr[iSoil], infiltParam.θ_Ini[iSoil]) # Not to have errors

					∑Infilt = INFILTRATION_MODEL(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt)

				elseif	option.infilt.OptimizeRun == "RunOptKs"  && option.θΨ ≠ "No"  #<>=<>=<>=<>=<>
					# Derive from laboratory
				
					SearchRange =[(log10(param.hydro.Ks_Min), log10(param.hydro.Ks_Max))]
					
					hydroInfilt.θr[iSoil] = min(hydroInfilt.θr[iSoil], infiltParam.θ_Ini[iSoil]) # Not to have errors

					Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(∑Infilt, ∑Infilt_Obs, hydroInfilt, infiltParam, iSoil, N_Infilt, T; Ks=10.0^P[1])[1]; SearchRange=SearchRange, NumDimensions=1, TraceMode=:silent)
		
						hydroInfilt.Ks[iSoil] = 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[1]

						∑Infilt = INFILTRATION_MODEL(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt)
				end # OptimizeRun = "Run"	
			end # iSoil

		# CONVERTING DIMENSIONS
			∑Infilt = CONVERT_INFILT_DIMENSIONS(hydroInfilt, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)

		return ∑Infilt, hydroInfilt
	end  # function: START_INFILTRATION

		
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_MODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_MODEL(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt)
			if option.infilt.Model == "Best_Univ" # <>=<>=<>=<>=<>
				bestUniv.BEST_UNIVERSAL_START(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt)		

			elseif option.infilt.Model == "QuasiExact" # <>=<>=<>=<>=<>
				# quasiExact.QUASIEXACT()

			end #  option.infilt.Model
			
			return ∑Infilt
		end  # function: INFILTRATION_MODEL


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERT_INFILT_DIMENSIONS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		function CONVERT_INFILT_DIMENSIONS(hydroInfilt, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)
			if option.infilt.Model == "Best_Univ" && option.infilt.SingleDoubleRing == "Double" && option.infilt.OutputDimension == "1D" # <>=<>=<>=<>=<>

				println("    ~ Converting $(option.infilt.Model) Infilt_3D => Infilt_1D ~")

				∑Infilt = bestUniv.CONVERT_3D_2_1D(hydroInfilt, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)

			elseif option.infilt.Model == "Best_Univ" && option.infilt.SingleDoubleRing == "Single" && option.infilt.OutputDimension == "3D" # <>=<>=<>=<>=<>

				println("    ~ Converting $(option.infilt.Model) Infilt_1D => Infilt_3D ~")

				∑Infilt = bestUniv.CONVERT_3D_2_1D(hydroInfilt, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)

			elseif option.infilt.Model == "QuasiExact" # <>=<>=<>=<>=<>
				# quasiExact.QUASIEXACT()
			end #  option.infilt

			return  ∑Infilt
		end  # function: CONVERT_INFILT_DIMENSIONS

	
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : HYDRO_PROCESS
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OBJECTIVE_FUNCTION(∑Infilt, ∑Infilt_Obs, hydroInfilt, infiltParam, iSoil, N_Infilt, T; σ=hydroInfilt.σ[iSoil], Ψm=hydroInfilt.Ψm[iSoil], θr=hydroInfilt.θr[iSoil], θs=hydroInfilt.θs[iSoil], Ks=hydroInfilt.Ks[iSoil])

				hydroInfilt.θs[iSoil] = θs
				hydroInfilt.θr[iSoil] = θr
				hydroInfilt.Ks[iSoil] = Ks
				hydroInfilt.σ[iSoil] = σ
				hydroInfilt.Ψm[iSoil] = Ψm

				hydroInfilt.θsMat[iSoil] = θs
				hydroInfilt.ΨmMac[iSoil] = Ψm
				hydroInfilt.σMac[iSoil] = σ


				∑Infilt = INFILTRATION_MODEL(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt)
				
				∑Of = 0
				for iT=1:N_Infilt[iSoil]
					∑Of += (∑Infilt[iSoil,iT] - ∑Infilt_Obs[iSoil,iT]) ^ 2.0					
				end

				return ∑Of
			end  # function: HYDRO_PROCESS
	

end  # module infiltration
