# =============================================================
#		MODULE: mainInfiltration
# =============================================================
module mainInfilt
	using ..option, ..sorptivity, ..best, ..param, ..wrc, ..kunsat, ...optInfilt
	export MAIN_INFILT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : MAIN_INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function MAIN_INFILT(N_SoilSelect, T, ∑Infilt, ∑Psd, N_Infilt, infilt, hydro)

		# COMPUTING INFILTRATION FROM PHYSICAL HYDRAULIC PARAMETERS
			if option.θΨ ≠ "No" && (option.infilt.OptimizeRun == "Run" ||  option.infilt.OptimizeRun == "RunOpt") && option.infilt.Model=="Simplified"
				global Infilt_Best_HydroObs = RUN_BEST(N_SoilSelect, infilt, hydro)	
			else
				global Infilt_Best_HydroObs = []
			end  # if: option.infilt.OptimizeRun == "Opt"

		# COMPUTING HYDRAULIC PARAMETERS FROM BEST
			if option.infilt.OptimizeRun  == "Opt" || option.infilt.OptimizeRun  == "RunOpt"		
				optInfilt.kg.OPT_INFILTRATION_BEST(N_SoilSelect, T, ∑Infilt, ∑Psd, N_Infilt, infilt)
			end

		return  Infilt_Best_HydroObs
	end  # function: MAIN_INFILT

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>		



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RUN_BEST
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RUN_BEST(N_SoilSelect, infilt, hydro)

			# We vary iSoil & SeIni & Number of time steps
			Infilt_Best_HydroObs = zeros(Float64, (N_SoilSelect, length(param.infilt.SeIni_Output), param.infilt.Npoint_Infilt))

			# For ever soil
			for iSoil=1:N_SoilSelect
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
			end  # for iSoil=1:N_SoilSelect
				
			return Infilt_Best_HydroObs
		end  # function: RUN_BEST



end  # module mainInfiltration
# ............................................................