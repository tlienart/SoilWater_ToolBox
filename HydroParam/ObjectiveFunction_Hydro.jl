module ofHydro 
	using ..option, ..stat, ..wrc, ..kunsat
	export  WRC_KUNSAT
	  
	function WRC_KUNSAT(iSoil, K_Kθ, N_KΨ, N_θΨ, θ_θΨ, Ψ_Kθ, Ψ_θΨ, hydro) 

		 # === OF θΨ ====
			θ_Obs = Array{Float64}(undef, N_θΨ[iSoil])
			θ_Sim = Array{Float64}(undef, N_θΨ[iSoil])
			
			@simd for iΨ in 1:N_θΨ[iSoil]
				θ_Obs[iΨ] = θ_θΨ[iSoil,iΨ]
				Ψ_Obs = Ψ_θΨ[iSoil,iΨ]
				θ_Sim[iΨ] = wrc.kg.Ψ_2_θDual(Ψ_Obs, hydro)
			end # for

			Of_θΨ = stat.NASH_SUTCLIFFE_ERRORmin(θ_Obs, θ_Sim; Power=2.0)


		 # === OF Kunsat ====
			Of_Kunsat = 0.0
			if option.KunsatΨ
				Kunsat_Obs_Ln = Array{Float64}(undef, N_KΨ[iSoil])
				Kunsat_Sim_Ln = Array{Float64}(undef, N_KΨ[iSoil])
				
				@simd for iΨ in 1:N_KΨ[iSoil]
					Kunsat_Obs_Ln[iΨ] = log1p(K_Kθ[iSoil,iΨ])
					Ψ_Obs =  Ψ_Kθ[iSoil,iΨ]
				
				Kunsat_Sim_Ln[iΨ] = log1p(kunsat.kg.Ψ_2_KUNSAT(Ψ_Obs, hydro))
				end # for

				Of_Kunsat = stat.NASH_SUTCLIFFE_ERRORmin(Kunsat_Obs_Ln, Kunsat_Sim_Ln)
			end #  option.KunsatΨ

		 return Of = 0.5 * Of_θΨ + 0.5 * Of_Kunsat

	end # WRC_KUNSAT

end # module ofHydaulic