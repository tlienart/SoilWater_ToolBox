# =============================================================
#		MODULE: mainInfiltration
# =============================================================
module mainInfilt
	using ..option, ..sorptivity, ..best, ..param
	export MAIN_INFILT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : MAIN_INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function MAIN_INFILT(N_SoilSelect, T, ∑Infilt, N_Infilt, infilt, hydro)


		# COMPUTING INFILTRATION FROM PHYSICAL HYDRAULIC PARAMETERS
		if option.θΨ ≠ "No" && (option.infilt.OptimizeRun == "Run" ||  option.infilt.OptimizeRun == "RunOpt")
			# Infilt: We vary iSoil & SeIni & Number of time steps
			Infilt_Best_HydroObs = zeros(Float64, (N_SoilSelect, length(param.infilt.SeIni_Output), 1000))
			N_ΔT_Best_HydroObs =  zeros(Int, (N_SoilSelect, length(param.infilt.SeIni_Output)))
			
			# For ever soil
			for iSoil=1:N_SoilSelect
				iSe_Ini = 0

				# For different initial Se_Ini defined by users param.infilt.SeIni_Output
				for Se_Ini in param.infilt.SeIni_Output
					iSe_Ini += 1

					# When to stop
						T_Max = best.TIME_TRANS_STEADY_INDEP(iSoil, Se_Ini, hydro, infilt) * param.infilt.Coeff_TransSteady

					# Number of time increment such that we stop one day!
						N_ΔT_Best_HydroObs[iSoil, iSe_Ini] = Int(T_Max ÷ param.infilt.ΔT_Output + 2)

						T = 0.0
						for iT=1:N_ΔT_Best_HydroObs[iSoil, iSe_Ini]
							T += param.infilt.ΔT_Output

							Infilt_Best_HydroObs[iSoil, iSe_Ini, iT] = best.BEST(iSoil, T, Se_Ini, hydro, infilt)
						end  # for iT=N_\Delta	
				end # for Se_Ini in param.infilt.SeIni_Output

			end  # for iSoil=1:N_SoilSelect
		else
			Infilt_Best_HydroObs = []
			N_ΔT_Best_HydroObs = []
		end  # if: option.infilt.OptimizeRun == "Opt"


		
		
		return Infilt_Best_HydroObs, N_ΔT_Best_HydroObs
	end  # function: MAIN_INFILT
	
end  # module mainInfiltration
# ............................................................