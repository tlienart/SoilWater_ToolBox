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
		function BOUNDARY_TOP_Ψ(discret, Flag_ReRun, hydro, iT, N_iZ, option, param, Q, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ)
			if Flag_ReRun
				Q[iT,2] = flux.Q!(option, discret, hydro, 2, iT, N_iZ, param, ΔHpond, ΔPr, ΔT, Ψ[iT,2], param.hyPix.Ψ_Top)
			else
				Q[iT,2] = flux.Q!(option, discret, hydro, 2, iT-1, N_iZ, param, ΔHpond, ΔPr, ΔT, Ψ[iT-1,2], param.hyPix.Ψ_Top)
			end
			
			θ[iT,1] = Ψ_2_θDual(option.hyPix, param.hyPix.Ψ_Top, 1, hydro)
			
			ΔPr[iT] = max((discret.ΔZ[1] * ((θ[iT,1] - θ[iT-1,1]) - hydro.So[1] * (param.hyPix.Ψ_Top - Ψ[iT-1,1]) * (θ[iT,1] / hydro.θs[1])) + ΔSink[iT-1,1]) / ΔT[iT] + Q[iT,2], 0.0)
			
		return ΔPr
		end  # function: BOUNDARY_TOP_\Psi
	# ------------------------------------------------------------------
	
end  # module: boundary
# ............................................................