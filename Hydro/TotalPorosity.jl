# =============================================================
#		MODULE: totalPorosity
# =============================================================
module Φ
	export  ρB_2_Φ

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Φ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρB_2_Φ(N_SoilSelect, RockW, ρ_Rock, ρbSoil, ρp_Fine)

			Φ = Array{Float64}(undef, (N_SoilSelect))

			for iSoil=1:N_SoilSelect
				# Φ_Fine = 1.0 - ρbSoil[iSoil] / ρp_Fine[iSoil]
				# [Volume Rock in soil] - [Fine earth bulk density]
				Φ_Fine = 1.0 - (RockW[iSoil] / ρ_Rock[iSoil]) - ((1.0 - RockW[iSoil]) *ρbSoil[iSoil]) / ((1.0 - RockW[iSoil] / ρ_Rock[iSoil]) * ρp_Fine[iSoil])


				# Vrock = RockW[iSoil] / ρ_Rock[iSoil]


				Φ[iSoil] = Φ_Fine
			end # for
			
			return Φ
		end  # function: Φ
	
end  # module: totalPorosity
# ............................................................