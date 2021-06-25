# =============================================================
#		MODULE: hydroRealation
# =============================================================
module hydroRelation
import BlackBoxOptim
import ..tool
export σ_2_Ψm⍰, σ_2_θr

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_Ψm⍰(iZ, hydro)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function σ_2_Ψm⍰(σ₁, Ψσ, Ψm_Min, Ψm_Max; Pσ=3.0)
			Ψm = Ψσ * exp(σ₁*Pσ)
			return Ψm = max(min(Ψm , Ψm_Max), Ψm_Min)
		end # function: σ_2_Ψm⍰


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_θr
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNCTION_σ_2_Ψm_SOFTWARE(hydro₂, iZ, option₂, param; Pσ=3.0)
			if (option₂.σ_2_Ψm⍰ == :Constrained)
				Ψm_Min = hydroRelation.σ_2_Ψm⍰(hydro₂.σ[iZ], param.hydro.kg.Ψσ_Min, hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ])

				Ψm_Max = hydroRelation.σ_2_Ψm⍰(hydro₂.σ[iZ], param.hydro.kg.Ψσ_Max, hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ])
				
				hydro₂.Ψm[iZ] = tool.norm.∇NORM_2_PARAMETER(hydro₂.Ψm[iZ], Ψm_Min, Ψm_Max)

			elseif (option₂.σ_2_Ψm⍰ == :UniqueRelationship) # <>=<>=<>=<>=<>
				hydro₂.Ψm[iZ] = hydroRelation.σ_2_Ψm⍰(hydro₂.σ[iZ], param.hydro.Ψσ, hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ])

			end #option.infilt.σ_2_Ψm⍰

		return hydro₂
		end


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_θr
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function σ_2_θr(hydro₂, iZ; α₁=17.5, α₂=4.0)
			σ_η = (hydro₂.σ[iZ] - hydro₂.σ_Min[iZ]) / (hydro₂.σ_Max[iZ] - hydro₂.σ_Min[iZ]) 
			
			return θr =  (hydro₂.θr_Max[iZ] * (1.0 - exp(-α₁ * σ_η ^ α₂))) / (1.0 - exp(-α₁ * 1.0 ^ α₂))
		end  # function: σ_2_θr
	
end  # module: hydroRealation
# ............................................................