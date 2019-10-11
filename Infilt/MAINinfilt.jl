# =============================================================
#		MODULE: mainInfiltration
# =============================================================
module mainInfilt
	using ..option, ..sorptivity, ..best, ..param, ..wrc, ..kunsat
	export MAIN_INFILT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : MAIN_INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function MAIN_INFILT(N_SoilSelect, T, ∑Infilt, N_Infilt, infilt, hydro)

		# For ever soil
		for iSoil=1:N_SoilSelect

			# COMPUTING INFILTRATION FROM PHYSICAL HYDRAULIC PARAMETERS
			if option.θΨ ≠ "No" && (option.infilt.OptimizeRun == "Run" ||  option.infilt.OptimizeRun == "RunOpt") && option.infilt.Model=="Simplified"
				global Infilt_Best_HydroObs = RUN_BEST(iSoil, N_SoilSelect, infilt, hydro, param)	
			else
				global Infilt_Best_HydroObs = []
			end  # if: option.infilt.OptimizeRun == "Opt"

		end  # for iSoil=1:N_SoilSelect
		
		
		return  Infilt_Best_HydroObs
	end  # function: MAIN_INFILT

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>		


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST_OPTIMISATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BEST_OPTIMISATION(N_SoilSelect, T, ∑Infilt, N_Infilt, infilt)
			
			return
		end  # function: BEST_OPTIMISATION



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RUN_BEST
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RUN_BEST(iSoil, N_SoilSelect, infilt, hydro, param)
			
			# We vary iSoil & SeIni & Number of time steps
			Infilt_Best_HydroObs = zeros(Float64, (N_SoilSelect, length(param.infilt.SeIni_Output), param.infilt.Npoint_Infilt))

			# For different initial Se_Ini defined by users param.infilt.SeIni_Output
			iSe_Ini = 0
			for iSe_Ini in 1: length(param.infilt.SeIni_Output)

				Se_Ini = param.infilt.SeIni_Output[iSe_Ini]
				# When to stop
					T_End = best.TIME_TRANS_STEADY_INDEP(iSoil, Se_Ini, hydro, infilt) * param.infilt.Coeff_TransSteady

					ΔT = (T_End + 2.0) / param.infilt.Npoint_Infilt # Time step

					T = 0.0
					for iT=1:param.infilt.Npoint_Infilt
						T += ΔT

						θ_Ini = wrc.Se_2_θ(Se_Ini, iSoil, hydro)

						Infilt_Best_HydroObs[iSoil, iSe_Ini, iT] = best.BEST(iSoil, T, θ_Ini, hydro, infilt)
					end  # for iT=N_\Delta	
			end # for Se_Ini in param.infilt.SeIni_Output
			
			return Infilt_Best_HydroObs
		end  # function: RUN_BEST



end  # module mainInfiltration
# ............................................................