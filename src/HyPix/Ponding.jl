
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
		function PONDING_SORPTIVITY!(discret, hydro, iT, option, param, Q, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ)
			
			Bparam = (2.0 - cst.β) / 3.0 + (Ψ_2_KUNSAT(option.hyPix, Ψ[iT-1,1], 1, hydro) / hydro.Ks[1]) * (1.0 + cst.β) / 3.0
			
			Infilt_Max =  (Sorptivity * √ΔT[iT] + Bparam * hydro.Ks[1] * ΔT[iT]) * param.hyPix.Cosα

			# Reduction of infiltration rate to avoid that too much water infiltrates into layer 1
				Infilt_Max = min(discret.ΔZ[1] * (hydro.θs[1] - θ[iT-1,1]) + ΔSink[iT,1], Infilt_Max)

			# ΔHpond is computed as
				ΔHpond[iT] = max(ΔPr[iT] + ΔHpond[iT-1] - Infilt_Max, 0.0)

				if isnan(ΔHpond[iT])
					error("ΔHpond = NaN")
				end
				
		return ΔHpond
		end  # function: PONDING
	#------------------------------------------------------------------
	
end  # module ponding
# ............................................................