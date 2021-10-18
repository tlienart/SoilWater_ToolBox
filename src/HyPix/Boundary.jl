# =============================================================
#		module: boundary
# =============================================================
module boundary
	import ..wrc: Ψ_2_θDual
	import ..flux

	export BOUNDARY_TOP_Ψ

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BOUNDARY_TOP_\Psi
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BOUNDARY_TOP_Ψ(discret, Flag_ReRun, hydro, iT, N_1, option, param, Q, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ)

			# If ponding ΔPr[iT] == 0.0
			if (ΔHpond[iT] == 0.0 && Flag_ReRun) || (ΔHpond[iT-1] == 0.0 && !Flag_ReRun)
				if !Flag_ReRun
					Q[iT,2] = flux.Q!(option, discret, hydro, 2, iT-1, N_1, param, ΔHpond, ΔPr, ΔT, θ, Ψ[iT-1,2], param.hyPix.Ψ_Top)
				else
					Q[iT,2] = flux.Q!(option, discret, hydro, 2, iT, N_1, param, ΔHpond, ΔPr, ΔT, θ, Ψ[iT,2], Ψ[iT,1])
				end
				
				θ[iT,1] = Ψ_2_θDual(option.hyPix, param.hyPix.Ψ_Top, 1, hydro)
				
				ΔPr[iT]  = (discret.ΔZ[1] * ((θ[iT,1] - θ[iT-1,1]) - hydro.So[1] * (param.hyPix.Ψ_Top- Ψ[iT-1,1]) * (θ[iT,1] / hydro.θs[1])) + ΔSink[iT,1]) / ΔT[iT] + Q[iT,2]

				ΔPr[iT]  =  (param.hyPix.Ψ_Top - Ψ[iT-1,1])

				ΔPr[iT] = max(ΔPr[iT], 0.0)
			else
				ΔPr[iT] == 0.0
			end

		return ΔPr
		end  # function: BOUNDARY_TOP_\Psi
	# ------------------------------------------------------------------
	
end  # module: boundary
# ............................................................