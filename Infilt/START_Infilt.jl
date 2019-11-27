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
			for iSoil=1:N_SoilSelect

				# No optimization required running from hydro
				if option.infilt.OptimizeRun == "Run" && option.θΨ ≠ "No"
					
					hydro.θr[iSoil] = min(hydro.θr[iSoil], infiltParam.θ_Ini[iSoil]) # Not to have errors
					∑Infilt = INFILTRATION_MODEL(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydro)
				end # OptimizeRun = "Run"	
			end

		# CONVERTING DIMENSIONS
			∑Infilt = CONVERT_INFILT_DIMENSIONS(hydro, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)

		return ∑Infilt
	end  # function: START_INFILTRATION


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

					println(∑Infilt)

				elseif option.infilt.Model == "Best_Univ" && option.infilt.SingleDoubleRing == "Single" && option.infilt.OutputDimension == "3D" # <>=<>=<>=<>=<>

					println("    ~ Converting $(option.infilt.Model) Infilt_1D => Infilt_3D ~")

					∑Infilt = bestUniv.CONVERT_3D_2_1D(hydroInfilt, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)


				elseif option.infilt.Model == "QuasiExact" # <>=<>=<>=<>=<>
					return
					# quasiExact.QUASIEXACT()

				end #  option.infilt
			end  # function: CONVERT_INFILT_DIMENSIONS

		# Sorptivity_Sum = 0.0
		# for iSoil=1:N_SoilSelect
		# 	

		# 	Sorptivity = sorptivity.SORPTIVITY(infiltParam.θ_Ini[iSoil], iSoil,  hydro)
		# 	println(iSoil," " ,Sorptivity)

		# 	Sorptivity_Sum += Sorptivity

		# end
		# println("Sum = $Sorptivity_Sum")

		# # COMPUTING INFILTRATION FROM PHYSICAL HYDRAULIC PARAMETERS
		# 	if option.θΨ ≠ "No" && (option.infilt.OptimizeRun == "Run" ||  option.infilt.OptimizeRun == "RunOpt") && option.infilt.Model=="Best_Univ" && option.infilt.Plot_SeIni_Range 

		# 		 ∑Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange = RUN_BEST_SeIni_RANGE(N_SoilSelect, infiltParam, hydro)	
		# 	else
		# 		 ∑Infilt_Best_HydroObs_SeIniRange = []
		# 		T_Best_HydroObs_SeIniRange = []
		# 	end  # if: option.infilt.OptimizeRun == "Opt"

		# 	if option.θΨ ≠ "No" && option.infilt.OptimizeRun  == "Run" || option.infilt.OptimizeRun  == "RunOpt" && option.infilt.Model=="Best_Univ"

		# 		∑Infilt_Best_HydroObs, Tinfilt_Best_HydroObs = RUN_BEST(N_SoilSelect, Tinfilt, N_Infilt, infiltParam, hydro)
		# 	else
		# 		 ∑Infilt_Best_HydroObs 	= []
		# 		Tinfilt_Best_HydroObs 	= []
		# 	end


		# # COMPUTING HYDRAULIC PARAMETERS FROM BEST
		# 	if option.infilt.OptimizeRun  == "Run" || option.infilt.OptimizeRun  == "RunOpt" && option.infilt.Plot_SeIni_Range	

		# 		optInfilt.kg.OPT_INFILTRATION_BEST(N_SoilSelect, Tinfilt, ∑Infilt, ∑Psd, N_Infilt, infiltParam)
		# 	end

		# return   ∑Infilt_Best_HydroObs, ∑Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange, Tinfilt_Best_HydroObs



	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>		

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RUN_BEST_SeIni_RANGE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RUN_BEST_SeIni_RANGE(N_SoilSelect, infiltParam, hydro)
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

						T_End = best.TIME_TRANS_STEADY_INDEP(iSoil, θ_Ini, hydro, infiltParam) * param.infilt.Coeff_TransSteady

						ΔT = (T_End + 2.0) / param.infilt.Npoint_Infilt # Time step

						# Cumulating time
						T_Best_HydroObs_SeIniRange[1] = ΔT
						for iT = 2:param.infilt.Npoint_Infilt
							T_Best_HydroObs_SeIniRange[iT] = T_Best_HydroObs_SeIniRange[iT-1] + ΔT
						end

					# looping
						for iT = 1:param.infilt.Npoint_Infilt
							θ_Ini = wrc.Se_2_θ(Se_Ini, iSoil, hydro)

							 ∑Infilt_Best_HydroObs_SeIniRange[iSoil, iSe_Ini, iT] = best.BEST(iSoil, "1D", Sorptivity, T_Best_HydroObs_SeIniRange[iT], θ_Ini, hydro, infiltParam)
						end  # for iT=N_\Delta	
				end # for Se_Ini in param.infilt.SeIni_Output	
			end  # for iSoil=1:N_SoilSelect
				
			return  ∑Infilt_Best_HydroObs_SeIniRange, T_Best_HydroObs_SeIniRange
		end  # function: RUN_BEST_SeIni_RANGE

end  # module infiltration
