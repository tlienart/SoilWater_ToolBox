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
				Φ_Fine = 1.0 - ρbSoil[iSoil] / ρp_Fine[iSoil]

				Φ[iSoil] = Φ_Fine
			end # for
			
			return Φ
		end  # function: Φ
	
end  # module: totalPorosity
# ............................................................