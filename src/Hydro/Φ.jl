# =============================================================
#		MODULE: totalPorosity
# =============================================================
module Φ
	export  ρB_2_Φ

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Φ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρB_2_Φ(N_SoilSelect, RockFragment, ρₚ_Rock, ρᵦ_Soil, ρₚ_Fine)

			Φ = Array{Float64}(undef, (N_SoilSelect))

			for iZ=1:N_SoilSelect
				# Vrock = RockFragment[iZ] / ρₚ_Rock[iZ]
				Φ[iZ] = 1.0 - (RockFragment[iZ] * ρᵦ_Soil[iZ] / ρₚ_Rock[iZ]) - ((1.0 - RockFragment[iZ]) *ρᵦ_Soil[iZ] / ρₚ_Fine[iZ])
			end # for
			
			return Φ
		end  # function: Φ
	
end  # module: totalPorosity
# ............................................................