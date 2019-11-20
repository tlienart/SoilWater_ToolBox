# =============================================================
#		MODULE: mainInfiltration
# =============================================================
module mainInfilt
	import ..option, ..sorptivity, ..best, ..param, ..wrc, ..kunsat, ...optInfilt
	export START_INFILT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function START_INFILT(N_SoilSelect, Tinfilt, ∑Infilt, ∑Psd, N_Infilt, infilt, hydro)

		# COMPUTING INFILTRATION FROM PHYSICAL HYDRAULIC PARAMETERS
			if option.θΨ ≠ "No" && (option.infilt.OptimizeRun == "Run" ||  option.infilt.OptimizeRun == "RunOpt") && option.infilt.Model=="Best_JJ" && option.infilt.SeIni_Range 

				 ∑Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange = RUN_BEST_SeIni_RANGE(N_SoilSelect, infilt, hydro)	
			else
				 ∑Infilt_Best_HydroObs_SeIniRange = []
				T_Best_HydroObs_SeIniRange = []
			end  # if: option.infilt.OptimizeRun == "Opt"

			if option.θΨ ≠ "No" && option.infilt.OptimizeRun  == "Run" || option.infilt.OptimizeRun  == "RunOpt" && option.infilt.Model=="Best_JJ"

				∑Infilt_Best_HydroObs, Tinfilt_Best_HydroObs = RUN_BEST(N_SoilSelect, Tinfilt, N_Infilt, infilt, hydro)
			else
				 ∑Infilt_Best_HydroObs 	= []
				Tinfilt_Best_HydroObs 	= []
			end


		# COMPUTING HYDRAULIC PARAMETERS FROM BEST
			if option.infilt.OptimizeRun  == "Run" || option.infilt.OptimizeRun  == "RunOpt" && option.infilt.SeIni_Range	

				optInfilt.kg.OPT_INFILTRATION_BEST(N_SoilSelect, Tinfilt, ∑Infilt, ∑Psd, N_Infilt, infilt)
			end

		return   ∑Infilt_Best_HydroObs, ∑Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange, Tinfilt_Best_HydroObs
	end  # function: START_INFILT


	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RUN_BEST
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RUN_BEST(N_SoilSelect, Tinfilt, N_Infilt, infilt, hydro)

			Tinfilt_Best_HydroObs = zeros(Float64, (N_SoilSelect, param.infilt.Npoint_Infilt))
			∑Infilt_Best_HydroObs = zeros(Float64, (N_SoilSelect, param.infilt.Npoint_Infilt) )
		
			# For ever soil
			for iSoil=1:N_SoilSelect
				# TIME CONTROL
					T_End = maximum(Tinfilt[iSoil,1:N_Infilt[iSoil]])
	
					ΔT = T_End / param.infilt.Npoint_Infilt # Time step

					# Cumulating time
					Tinfilt_Best_HydroObs[1] = ΔT
					for iT = 2:param.infilt.Npoint_Infilt
						Tinfilt_Best_HydroObs[iSoil,iT] = Tinfilt_Best_HydroObs[iSoil,iT-1] + ΔT
					end

				# SORPTIVITY
					Sorptivity = sorptivity.SORPTIVITY(infilt.θ_Ini[iSoil], iSoil, hydro)
				
				# looping
					for iT = 1:param.infilt.Npoint_Infilt
						 ∑Infilt_Best_HydroObs[iSoil,iT] = best.BEST(iSoil, option.infilt.Dimension, Sorptivity, Tinfilt_Best_HydroObs[iSoil,iT], infilt.θ_Ini[iSoil], hydro, infilt)
					end  # for iT=N_\Delta	
			end  # for iSoil=1:N_SoilSelect
				
			return  ∑Infilt_Best_HydroObs, Tinfilt_Best_HydroObs
		end # function: RUN_BEST


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RUN_BEST_SeIni_RANGE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RUN_BEST_SeIni_RANGE(N_SoilSelect, infilt, hydro)
			# We vary iSoil & SeIni & Number of time steps
			 ∑Infilt_Best_HydroObs_SeIniRange = zeros(Float64, (N_SoilSelect, length(param.infilt.SeIni_Output), param.infilt.Npoint_Infilt))
			
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

							 ∑Infilt_Best_HydroObs_SeIniRange[iSoil, iSe_Ini, iT] = best.BEST(iSoil, "1D", Sorptivity, T_Best_HydroObs_SeIniRange[iT], θ_Ini, hydro, infilt)
						end  # for iT=N_\Delta	
				end # for Se_Ini in param.infilt.SeIni_Output	
			end  # for iSoil=1:N_SoilSelect
				
			return  ∑Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange
		end  # function: RUN_BEST_SeIni_RANGE

end  # module mainInfiltration
# ............................................................