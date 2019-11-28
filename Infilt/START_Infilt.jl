# =============================================================
#		MODULE: infiltration
# =============================================================
module infilt
	import ..option, ..sorptivity, ..param, ..wrc, ..kunsat, ...opt, ..infiltInitialize, ..bestUniv
	export START_INFILTRATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function START_INFILTRATION(N_SoilSelect, Tinfilt, ∑Infilt_Obs, ∑Psd, N_Infilt, infiltParam, hydro)

		# INITIALIZE
			∑Infilt, T, hydroInfilt = infiltInitialize.INFILT_INITIALIZE(N_SoilSelect, ∑Psd, infiltParam, Tinfilt, N_Infilt)

		# OPTIONS
			# For every soils
			for iSoil=1:N_SoilSelect

				# No optimization required running from hydro
				if option.infilt.OptimizeRun == "Run" && option.θΨ ≠ "No"
					hydroInfilt = hydro

					hydroInfilt.θr[iSoil] = min(hydroInfilt.θr[iSoil], infiltParam.θ_Ini[iSoil]) # Not to have errors

					∑Infilt = INFILTRATION_MODEL(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt)

				elseif	option.infilt.OptimizeRun == "RunOptKs"  && option.θΨ ≠ "No"
					hydroInfilt = hydro

					hydroInfilt.θr[iSoil] = min(hydroInfilt.θr[iSoil], infiltParam.θ_Ini[iSoil]) # Not to have errors

					∑Infilt = INFILTRATION_MODEL(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt)
				end # OptimizeRun = "Run"	
			end

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

	

	

end  # module infiltration
