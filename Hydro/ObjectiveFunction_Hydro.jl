module ofHydro 
	import ..option, ...stats, ..wrc, ..kunsat
	export  WRC_KUNSAT
	  
	function OF_WRC_KUNSAT(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; Wof = 0.5) 

		 # === OF θΨ ====
			θ_Obs = Array{Float64}(undef, N_θΨ[iSoil])
			θ_Sim = Array{Float64}(undef, N_θΨ[iSoil])
			
			for iΨ = 1:N_θΨ[iSoil]
				θ_Obs[iΨ] = θ_θΨ[iSoil,iΨ]
				θ_Sim[iΨ] = wrc.Ψ_2_θDual(Ψ_θΨ[iSoil,iΨ], iSoil, hydro)
			end # for iΨ in 1:N_θΨ[iSoil]

			Of_θΨ = stats.NASH_SUTCLIFE_MINIMIZE(θ_Obs, θ_Sim)


		 # === OF Kunsat ====
			if option.KunsatΨ
				Kunsat_Obs_Ln = Array{Float64}(undef, N_KΨ[iSoil])
				Kunsat_Sim_Ln = Array{Float64}(undef, N_KΨ[iSoil])
				for iΨ = 1:N_KΨ[iSoil]
					Kunsat_Obs_Ln[iΨ] = log1p(K_KΨ[iSoil,iΨ])
						
					Kunsat_Sim_Ln[iΨ] = log1p(kunsat.Ψ_2_KUNSAT(Ψ_KΨ[iSoil,iΨ], iSoil, hydro))
				end # for iΨ in 1:N_KΨ[iSoil]

				Of_Kunsat = stats.NASH_SUTCLIFE_MINIMIZE(Kunsat_Obs_Ln, Kunsat_Sim_Ln)
			end #  option.KunsatΨ

			Of = Wof * Of_θΨ + (1.0 - Wof) * Of_Kunsat

		 return Of, Of_θΨ, Of_Kunsat

	end # function OF_WRC_KUNSAT

end # module ofHydaulic