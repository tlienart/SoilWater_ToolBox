module ofHydrolab
	import..stats, ..wrc, ..kunsat
	export  OF_WRC_KUNSAT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_WRC_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
		function OF_WRC_KUNSAT(optionₘ, iZ, θ_θΨobs, Ψ_θΨobs, N_θΨobs, K_KΨobs, Ψ_KΨobs, N_KΨobs, hydro, optim; Wof=0.5) 

			# === OF θΨ ====
				θ_Obs = fill(0.0::Float64, N_θΨobs[iZ])
				θ_Sim = fill(0.0::Float64, N_θΨobs[iZ])

				for iΨ = 1:N_θΨobs[iZ]
					θ_Obs[iΨ] = θ_θΨobs[iZ,iΨ]
					θ_Sim[iΨ] = wrc.Ψ_2_θDual(optionₘ,  Ψ_θΨobs[iZ,iΨ], iZ, hydro)
				end # for iΨ = 1:N_θΨobs[iZ]

				Of_θΨ = stats.NSE_MINIMIZE(θ_Obs[1:N_θΨobs[iZ]], θ_Sim[1:N_θΨobs[iZ]])

			# === OF Kunsat ====
			if "Ks" ∈ optim.ParamOpt
				if "Ks" ∈ optim.ParamOpt
					iStart = 1
				else
					iStart = 2
				end

				# Of_Kunsat = 0.0
				Kunsat_Obs_Ln = fill(0.0::Float64, N_KΨobs[iZ])
				Kunsat_Sim_Ln = fill(0.0::Float64, N_KΨobs[iZ])

				for iΨ = iStart:N_KΨobs[iZ]
					Kunsat_Obs_Ln[iΨ] = log1p(K_KΨobs[iZ,iΨ])
						
					Kunsat_Sim_Ln[iΨ] = log1p(kunsat.Ψ_2_KUNSAT(optionₘ, Ψ_KΨobs[iZ,iΨ], iZ, hydro))
				end # for iΨ = 1:N_KΨobs[iZ]

				Of_Kunsat = stats.NSE_MINIMIZE(Kunsat_Obs_Ln[iStart:N_KΨobs[iZ]], Kunsat_Sim_Ln[iStart:N_KΨobs[iZ]])			
				Of = Wof * Of_θΨ + (1.0 - Wof) * Of_Kunsat

			else		
				Of = Of_θΨ
				Of_Kunsat = 0.0
			end #  "Ks" ∈ optim.ParamOpt
		return Of, Of_θΨ, Of_Kunsat
		end # function OF_WRC_KUNSAT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_RMSE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
		function OF_RMSE(option, optionₘ, iZ, θ_θΨobs, Ψ_θΨobs, N_θΨobs, K_KΨobs, Ψ_KΨobs, N_KΨobs, hydro, optim) 

		# === OF θΨ ====
			θ_Obs = fill(0.0::Float64, N_θΨobs[iZ])
			θ_Sim = fill(0.0::Float64, N_θΨobs[iZ])

			for iΨ = 1:N_θΨobs[iZ]
				θ_Obs[iΨ] = θ_θΨobs[iZ,iΨ]
				θ_Sim[iΨ] = wrc.Ψ_2_θDual(optionₘ, Ψ_θΨobs[iZ,iΨ], iZ, hydro)
			end # for iΨ = 1:N_θΨobs[iZ]

			Rmse_θΨ = stats.RMSE(θ_Obs[1:N_θΨobs[iZ]], θ_Sim[1:N_θΨobs[iZ]])

		# === OF Kunsat ====
			if "Ks" ∈ optim.ParamOpt ||option.run.HydroLabθΨ⍰=="Run"
				if  "Ks" ∈ optim.ParamOpt
					iStart = 1
				else
					iStart = 2
				end

				Kunsat_Obs_Ln = fill(0.0::Float64, N_KΨobs[iZ])
				Kunsat_Sim_Ln = fill(0.0::Float64, N_KΨobs[iZ])

				for iΨ = iStart:N_KΨobs[iZ]
					Kunsat_Obs_Ln[iΨ] = log1p(K_KΨobs[iZ,iΨ])
						
					Kunsat_Sim_Ln[iΨ] = log1p(kunsat.Ψ_2_KUNSAT(optionₘ, Ψ_KΨobs[iZ,iΨ], iZ, hydro))
				end # for iΨ = 1:N_KΨobs[iZ]

				Rmse_KΨ = stats.RMSE(Kunsat_Obs_Ln[iStart:N_KΨobs[iZ]], Kunsat_Sim_Ln[iStart:N_KΨobs[iZ]])

				Rmse = (Rmse_θΨ + Rmse_KΨ) * 0.5
			else		
				Rmse = Rmse_θΨ
				Rmse_KΨ = 0.0
			end #  "Ks" ∈ optim.ParamOpt

	return Rmse, Rmse_KΨ, Rmse_θΨ
	end # OF_RMSE

end # module ofHydaulic
# ............................................................