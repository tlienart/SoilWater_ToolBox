# =============================================================
#		MODULE: totalPorosity
# =============================================================
module Φ
	export  ρᵦ_2_Φ

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Φ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρᵦ_2_Φ(N_SoilSelect, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)
			Φ = fill(0.0::Float64, N_SoilSelect)

			for iZ=1:N_SoilSelect
				if option.rockFragment.RockInjectedIncluded == :Included 
					Φ[iZ] = 1.0 - (RockFragment[iZ] * ρᵦ_Soil[iZ] / ρₚ_Rock[iZ]) - ((1.0 - RockFragment[iZ]) * ρᵦ_Soil[iZ] / ρₚ_Fine[iZ])

				elseif option.rockFragment.RockInjectedIncluded == :Injected 
					Φ[iZ] = (1.0 - ρᵦ_Soil[iZ] / ρₚ_Rock[iZ]) * (1.0 - RockFragment[iZ])

				end
			end # for
			
		return Φ
		end  # function: Φ
	
end  # module: totalPorosity
# ............................................................