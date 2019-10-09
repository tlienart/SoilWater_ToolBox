# =============================================================
#		MODULE: mainInfiltration
# =============================================================
module mainInfilt
	using ..option, ..sorptivity
	export MAIN_INFILT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : MAIN_INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function MAIN_INFILT(N_SoilSelect, T, ∑Infilt, N_Infilt, infilt, hydro)

		if option.θΨ ≠ "No" && (option.infilt.OptimizeRun == "Run" ||  option.infilt.OptimizeRun == "RunOpt")
			
			for iSoil=1:N_SoilSelect

				Sorptivity = sorptivity.kg.SORPTIVITY(iSoil, infilt.Se_Ini[iSoil], hydro)

				println(Sorptivity)
					
			end  # for iSoil=1:N_SoilSelect
			
			
		end  # if: option.infilt.OptimizeRun == "Opt"
		
		return
	end  # function: MAIN_INFILT
	
end  # module mainInfiltration
# ............................................................