
# =============================================================
#		MODULE: ponding
# =============================================================
module ponding
	import  ..cst
	import ..kunsat:Ψ_2_KUNSAT
	export PONDING_RUNOFF_SORPTIVITY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PONDING_RUNOFF_SORPTIVITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PONDING_RUNOFF_SORPTIVITY(discret, Hpond, hydro, iT::Int64, optionHypix, paramHypix, Sorptivity, ΔPr, ΔRunoff, ΔSink, ΔT, θ, Ψ)
			
			# PONDING
				Bparam = (2.0 - cst.β) / 3.0 + (Ψ_2_KUNSAT(optionHypix, Ψ[iT-1,1], 1, hydro) / hydro.Ks[1]) * (1.0 + cst.β) / 3.0
				
				Infilt_Max =  (Sorptivity * √ΔT[iT] + Bparam * hydro.Ks[1] * ΔT[iT]) * paramHypix.Cosα

				# Reduction of infiltration rate to avoid that too much water infiltrates into layer 1
					Infilt_Max = min(discret.ΔZ[1] * (hydro.θs[1] - θ[iT-1,1]) + ΔSink[iT,1], Infilt_Max)

					Hpond[iT] = max(ΔPr[iT] + Hpond[iT-1] - Infilt_Max, 0.0::Float64)
							
			# RUNOFF
				Hpond_Max = paramHypix.Cosα * paramHypix.Hpond_Max
				if Hpond[iT] > Hpond_Max
               ΔRunoff[iT] = Hpond[iT] - Hpond_Max
               Hpond[iT]   = Hpond_Max
				else
					ΔRunoff[iT] = 0.0::Float64
				end

		return Hpond, ΔRunoff
		end  # function: PONDING
	#------------------------------------------------------------------
	
end  # module ponding
# ............................................................