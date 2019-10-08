module ofHydro 
	using ..option, ...stats, ..wrc, ..kunsat
	export  WRC_KUNSAT
	  
	function OF_WRC_KUNSAT(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro) 

		 # === OF θΨ ====
			θ_Obs = Array{Float64}(undef, N_θΨ[iSoil])
			θ_Sim = Array{Float64}(undef, N_θΨ[iSoil])
			
			@simd for iΨ in 1:N_θΨ[iSoil]
				θ_Obs[iΨ] = θ_θΨ[iSoil,iΨ]
				Ψ_Obs = Ψ_θΨ[iSoil,iΨ]
				θ_Sim[iΨ] = wrc.kg.Ψ_2_θDual(Ψ_Obs, iSoil, hydro)
			end # for iΨ in 1:N_θΨ[iSoil]

			Of_θΨ = stats.NASH_SUTCLIFFE_ERRORmin(θ_Obs, θ_Sim; Power=2.0)


		 # === OF Kunsat ====
			Of_Kunsat = 0.0
			if option.KunsatΨ
				Kunsat_Obs_Ln = Array{Float64}(undef, N_KΨ[iSoil])
				Kunsat_Sim_Ln = Array{Float64}(undef, N_KΨ[iSoil])
				
				@simd for iΨ in 1:N_KΨ[iSoil]
					Kunsat_Obs_Ln[iΨ] = log1p(K_KΨ[iSoil,iΨ])
					Ψ_Obs =  Ψ_KΨ[iSoil,iΨ]
				
					Kunsat_Sim_Ln[iΨ] = log1p(kunsat.kg.Ψ_2_KUNSAT(Ψ_Obs, iSoil, hydro))
				end # for iΨ in 1:N_KΨ[iSoil]

				Of_Kunsat = stats.NASH_SUTCLIFFE_ERRORmin(Kunsat_Obs_Ln, Kunsat_Sim_Ln)
			end #  option.KunsatΨ

			Of = 0.5 * Of_θΨ + 0.5 * Of_Kunsat

		 return Of, Of_θΨ, Of_Kunsat

	end # function OF_WRC_KUNSAT

end # module ofHydaulic