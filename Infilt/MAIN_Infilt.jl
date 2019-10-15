# =============================================================
#		MODULE: mainInfiltration
# =============================================================
module mainInfilt
	using ..option, ..sorptivity, ..best, ..param, ..wrc, ..kunsat, ...optInfilt
	export MAIN_INFILT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : MAIN_INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function MAIN_INFILT(N_SoilSelect, Tinfilt, ∑Infilt, ∑Psd, N_Infilt, infilt, hydro)

		# COMPUTING INFILTRATION FROM PHYSICAL HYDRAULIC PARAMETERS
			if option.θΨ ≠ "No" && (option.infiltration.OptimizeRun == "Run" ||  option.infiltration.OptimizeRun == "RunOpt") && option.infiltration.Model=="Simplified"

				Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange = RUN_BEST_SeIni_RANGE(N_SoilSelect, infilt, hydro)	
			else
				Infilt_Best_HydroObs_SeIniRange = []
				T_Best_HydroObs_SeIniRange = []
			end  # if: option.infiltration.OptimizeRun == "Opt"

			# if option.infiltration.OptimizeRun  == "Run" || option.infiltration.OptimizeRun  == "RunOpt" 
			# 	RUN_BEST(N_SoilSelect, Tinfilt, ∑Infilt, N_Infilt, infilt, hydro)
			# end


		# COMPUTING HYDRAULIC PARAMETERS FROM BEST
			if option.infiltration.OptimizeRun  == "Run" || option.infiltration.OptimizeRun  == "RunOpt" && option.infiltration.SeIni_Range		
				optInfilt.kg.OPT_INFILTRATION_BEST(N_SoilSelect, Tinfilt, ∑Infilt, ∑Psd, N_Infilt, infilt)
			end

		return  Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange
	end  # function: MAIN_INFILT

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RUN_BEST
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RUN_BEST(N_SoilSelect, Tinfilt, ∑Infilt, N_Infilt, infilt, hydro)

			# Infilt_Best_HydroOb = zeros(Float64, (N_SoilSelect, param.infilt.Npoint_Infilt))
			Tinfilt_Best_HydroObs = zeros(Float64, (param.infilt.Npoint_Infilt))
		
			# For ever soil
			for iSoil=1:N_SoilSelect
				# For different initial Se_Ini defined by users param.infilt.SeIni_Output


				# TIME CONTROL
					T_End = maximum(Tinfilt)

					ΔT = (T_End + 2.0) / param.infilt.Npoint_Infilt # Time step

					# Cumulating time
					Tinfilt_Best_HydroObs[1] = ΔT
					for iT = 2:param.infilt.Npoint_Infilt
						Tinfilt_Best_HydroObs[iT] = Tinfilt_Best_HydroObs[iT-1] + ΔT
					end

				# SORPTIVITY
					Sorptivity = sorptivity.SORPTIVITY(infilt.θ_Ini[iSoil], iSoil, hydro)
				
				# looping
					for iT = 1:param.infilt.Npoint_Infilt

						Infilt_Best_HydroObs[iSoil,iT] = best.BEST(iSoil, option.infiltration.Dimension, Sorptivity, Tinfilt_Best_HydroObs[iT], infilt.θ_Ini[iSoil], hydro, infilt)
					end  # for iT=N_\Delta	
			end  # for iSoil=1:N_SoilSelect
				
			return Infilt_Best_HydroObs
		end # function: RUN_BEST


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RUN_BEST_SeIni_RANGE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RUN_BEST_SeIni_RANGE(N_SoilSelect, infilt, hydro)
			# We vary iSoil & SeIni & Number of time steps
			Infilt_Best_HydroObs_SeIniRange = zeros(Float64, (N_SoilSelect, length(param.infilt.SeIni_Output), param.infilt.Npoint_Infilt))
			
			# For ever soil
			for iSoil=1:N_SoilSelect
				# For different initial Se_Ini defined by users param.infilt.SeIni_Output
				iSe_Ini = 0
				for iSe_Ini in 1: length(param.infilt.SeIni_Output)

					Se_Ini = param.infilt.SeIni_Output[iSe_Ini]

					θ_Ini = wrc.Se_2_θ(Se_Ini, iSoil, hydro)

					Sorptivity = sorptivity.SORPTIVITY(θ_Ini, iSoil, hydro)
					
					# TIME CONTROL
						global T_Best_HydroObs_SeIniRange = zeros(Float64, param.infilt.Npoint_Infilt)

						T_End = best.TIME_TRANS_STEADY_INDEP(iSoil, θ_Ini, hydro, infilt) * param.infilt.Coeff_TransSteady

						ΔT = (T_End + 2.0) / param.infilt.Npoint_Infilt # Time step

						# Cumulating time
						T_Best_HydroObs_SeIniRange[1] = ΔT
						for iT = 2:param.infilt.Npoint_Infilt
							T_Best_HydroObs_SeIniRange[iT] = T_Best_HydroObs_SeIniRange[iT-1] + ΔT
						end

					# looping
						for iT = 1:param.infilt.Npoint_Infilt
							θ_Ini = wrc.Se_2_θ(Se_Ini, iSoil, hydro)

							Infilt_Best_HydroObs_SeIniRange[iSoil, iSe_Ini, iT] = best.BEST(iSoil, "1D", Sorptivity, T_Best_HydroObs_SeIniRange[iT], θ_Ini, hydro, infilt)
						end  # for iT=N_\Delta	
				end # for Se_Ini in param.infilt.SeIni_Output	
			end  # for iSoil=1:N_SoilSelect
				
			return Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange
		end  # function: RUN_BEST_SeIni_RANGE

end  # module mainInfiltration
# ............................................................