module stats
	import ..wrc, ..param
	export NASH_SUTCLIFE_MINIMIZE, NASH_SUTCLIFFE_θΨ, RELATIVE_ERR, NASH_SUTCLIFFE_EFFICIENCY
	using Statistics

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NASH_SUTCLIFE_MINIMIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NASH_SUTCLIFE_MINIMIZE(Obs, Sim; Power=2.0)
			N = length(Obs)
			Obs_Mean = Statistics.mean(Obs[1:N])
			Obs_Mean_Err = Statistics.sum(abs.(Obs_Mean .- Obs[1:N]) .^ Power)

			if Obs_Mean_Err < 0.000000000000001   
				Obs_Mean_Err = 1.0
			end
			
			Err = 0.0
			for i = 1:N
				Err += sum(abs(Sim[i] - Obs[i]) ^ Power)
			end
			return Err / Obs_Mean_Err 
		end  # function: NASH_SUTCLIFE_MINIMIZE

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NASH_SUTCLIFFE_θΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Data, Ψ_Sim, θ_Sim, hydro)
			Nse = zeros(Float64, N_SoilSelect)

			for iSoil = 1:N_SoilSelect	
				θΨ = zeros(Float64, N_Data[iSoil])
				for iRpart = 1:N_Data[iSoil]
					θΨ[iRpart] = wrc.Ψ_2_θDual(Ψ_Sim[iSoil,iRpart], iSoil, hydro)
				end
				Nse[iSoil] = 1.0 - NASH_SUTCLIFE_MINIMIZE(θΨ[1:N_Data[iSoil]], θ_Sim[iSoil,1:N_Data[iSoil]])	
			end

			# Cumulating the objective function to get the overview
			Nse_Mean = round(Statistics.mean(max.(Nse[1:N_SoilSelect],0.0)), digits=3)  # in case of negative value then it is set to 0
			Nse_Std  = round(Statistics.std(max.(Nse[1:N_SoilSelect],0.0)), digits=3)   # in case of negative value then it is set to 0

			return Nse, Nse_Mean, Nse_Std
		end # function NASH_SUTCLIFFE_θΨ

		
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NASH_SUTCLIFFE_EFFICIENCY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NASH_SUTCLIFFE_EFFICIENCY(;Obs=Obs, Sim=Sim, Power=2.0)
			return Nse = 1 - NASH_SUTCLIFE_MINIMIZE(Obs, Sim; Power=Power)
		end  # function: NASH_SUTCLIFFE_EFFICIENCY


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RELATIVE_ERR
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RELATIVE_ERR(;Obs=Obs, Sim=Sim)
			return Err = 1.0 - abs(Obs - Sim) / Obs
		end

end # module