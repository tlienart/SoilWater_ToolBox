
# =============================================================
#		MODULE: ponding
# =============================================================
module ponding
	import  ..cst
	import ..kunsat:Ψ_2_KUNSAT
	export PONDING_SORPTIVITY!

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PONDING_SORPTIVITY!
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PONDING_SORPTIVITY!(discret, hydro, iT, optionHypix, paramHypix, Q, Sorptivity, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ)
			
			Bparam = (2.0 - cst.β) / 3.0 + (Ψ_2_KUNSAT(optionHypix, Ψ[iT-1,1], 1, hydro) / hydro.Ks[1]) * (1.0 + cst.β) / 3.0
			
			Infilt_Max =  (Sorptivity * √ΔT[iT] + Bparam * hydro.Ks[1] * ΔT[iT]) * paramHypix.Cosα

			# Reduction of infiltration rate to avoid that too much water infiltrates into layer 1
				Infilt_Max = min(discret.ΔZ[1] * (hydro.θs[1] - θ[iT-1,1]) + ΔSink[iT,1], Infilt_Max)

			# Hpond is computed as
				Hpond[iT] = max(ΔPr[iT] + Hpond[iT-1] - Infilt_Max, 0.0)

				if isnan(Hpond[iT])
					error("Hpond = NaN")
				end
				
		return Hpond
		end  # function: PONDING
	#------------------------------------------------------------------
	
end  # module ponding
# ............................................................