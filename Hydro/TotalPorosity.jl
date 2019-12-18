# =============================================================
#		MODULE: totalPorosity
# =============================================================
module Φ
	export  ρB_2_Φ

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Φ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρB_2_Φ(N_SoilSelect, ρb_Fine, ρp_Fine, ρb_Rock, ρp_Rock)

			Φ = Array{Float64}(undef, (N_SoilSelect))

			for iSoil=1:N_SoilSelect
				Φ_Fine = 1.0 - ρb_Fine[iSoil] / ρp_Fine[iSoil]

				if ρb_Rock[iSoil] > 0.0 && ρp_Rock[iSoil] > 0.0
					Φ_Rock = 1.0 - ρb_Rock[iSoil] / ρp_Rock[iSoil]
				else
					Φ_Rock = 0.0
				end

				Φ[iSoil] = Φ_Fine

			end # for
			
			return Φ
		end  # function: Φ
	
end  # module: totalPorosity
# ............................................................