module stats
	import ..wrc, ..kunsat
	export RMSE, NSE_MINIMIZE, NSE_θΨ, RELATIVE_ERR, NSE_EFFICIENCY, LINEAR_REGRESSION, NSE_WILMOT, RMSE_CONCORDANCE_CORELATION_COEFICIENT
	using Statistics

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NSE_MINIMIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RMSE(Obs, Sim; Power=2.0)
			N = length(Obs)
			Err = 0.0
			iCount = 1
			
			for i = 1:N
				if !(isnan(Sim[i]) || isnan(Obs[i]))
					Err += abs(Sim[i] - Obs[i]) ^ Power
					iCount += 1
				end
			end
		return (Err / Float64(iCount)) ^ (1.0 / Power)
		end # function RMSE


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NSE_MINIMIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NSE(Obs, Sim; Power=2.0)
			N = length(Obs)
			Err = 0.0
			iCount = 1

			∑Obs = 0.0
			iCount = 1
			for i = 1:N
				if !(isnan(Sim[i]) || isnan(Obs[i]))
					∑Obs += Obs[i]
					iCount += 1
				end
			end
			Obs_Mean = ∑Obs / Float64(iCount)
			
			Obs_Mean_Err = 0.0
			for i = 1:N
				if !(isnan(Sim[i]) || isnan(Obs[i]))
					Err += abs(Sim[i] - Obs[i]) ^ Power

					Obs_Mean_Err += abs(Obs_Mean - Obs[i]) ^ Power
				end
			end
		return 1.0 - (Err / Obs_Mean_Err)
		end # function NSE


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : NSE_WILMOT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function NSE_WILMOT(Obs, Sim; C=2)
				N = length(Obs)
				∑Obs = 0.0
				iCount = 1
			
				for i = 1:N
					if !(isnan(Sim[i]) || isnan(Obs[i]))
						∑Obs += Obs[i]
						iCount += 1
					end
				end
				Obs_Mean = ∑Obs / Float64(iCount)
				
				# or can we used just this -->   Obs_Mean = Statistics.mean(Obs[1:N])
				Nse_Wilmot = 0.0
				for i = 1:N
					if Statistics.sum( abs.(Sim[1:N] .- Obs[1:N]) ) ≤ C .* Statistics.sum( abs.(Obs[1:N] .- Obs_Mean) )
						Nse_Wilmot = 1.0 .- (Statistics.sum( abs.(Sim[1:N] .- Obs[1:N]))) / (C * (Statistics.sum(abs.(Obs[1:N] .- Obs_Mean))))
					else
						Nse_Wilmot = (C * (Statistics.sum(abs.(Obs[1:N] .- Obs_Mean)))) / (Statistics.sum(abs.(Sim[1:N] .- Obs[1:N]))) - 1.0 
					end
				end # for

			return Nse_Wilmot
			end  # function: NSE_WILMOT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NSE_CONCORDANCE_CORELATION_COEFICIENT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NSE_CONCORDANCE_CORELATION_COEFICIENT(Obs, Sim)
			PersonCorelation = Statistics.cor(Obs, Sim)
			σₒᵦₛ             = Statistics.std(Obs)
			σₛᵢₘ             = Statistics.std(Sim)
			Meanₒᵦₛ          = Statistics.mean(Obs)
			Meanₛᵢₘ          = Statistics.mean(Sim)

		return  (2.0 * PersonCorelation * σₒᵦₛ * σₛᵢₘ) / (σₒᵦₛ ^ 2.0 +  σₛᵢₘ ^ 2.0 + (Meanₒᵦₛ - Meanₛᵢₘ) ^ 2.0)
		end  # function: BestConcordanceCorrelationCoeficient
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NSE_CONCORDANCE_CORELATION_COEFICIENT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RMSE_CONCORDANCE_CORELATION_COEFICIENT(Obs, Sim)
			PersonCorelation = Statistics.cor(Obs, Sim)
			σₒᵦₛ             = Statistics.std(Obs)
			σₛᵢₘ             = Statistics.std(Sim)
			Meanₒᵦₛ          = Statistics.mean(Obs)
			Meanₛᵢₘ          = Statistics.mean(Sim)

		Ccc = 1.0 - (2.0 * PersonCorelation * σₒᵦₛ * σₛᵢₘ) / (σₒᵦₛ ^ 2.0 +  σₛᵢₘ ^ 2.0 + (Meanₒᵦₛ - Meanₛᵢₘ) ^ 2.0)

		return max(Ccc, 0.0)
		end  # function: BestConcordanceCorrelationCoeficient
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NSE_MINIMIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NSE_MINIMIZE(Obs, Sim; Power=2.0)
			N = length(Obs)
			Obs_Mean = Statistics.mean(Obs[1:N])
			Obs_Mean_Err = Statistics.sum(abs.(Obs_Mean .- Obs[1:N]) .^ Power)

			if Obs_Mean_Err < 0.000000000000001   
				Obs_Mean_Err = 1.0
			end
			
			Err = 0.0
			for i = 1:N
				Err += abs(Sim[i] - Obs[i]) ^ Power
			end
		return Err / Obs_Mean_Err 
		end  # function: NSE_MINIMIZE

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NSE_θΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NSE_θΨ(hydro, N_Data, NiZ, optionₘ, θobs, Ψobs)
			Nse = zeros(Float64, NiZ)

			for iZ = 1:NiZ	
				θΨsim = zeros(Float64, N_Data[iZ])
				for iΨ = 1:N_Data[iZ]
					θΨsim[iΨ] = wrc.Ψ_2_θDual(optionₘ, Ψobs[iZ,iΨ], iZ, hydro)
				end
				Nse[iZ] = 1.0 - NSE_MINIMIZE( θobs[iZ,1:N_Data[iZ]], θΨsim[1:N_Data[iZ]])	
			end

			# Cumulating the objective function to get the overview
			Nse_Mean = round(Statistics.mean(max.(Nse[1:NiZ] , 0.0)), digits=3)  # in case of negative value then it is set to 0
			Nse_Std  = round(Statistics.std(max.(Nse[1:NiZ] , 0.0)), digits=3)   # in case of negative value then it is set to 0

		return Nse, Nse_Mean, Nse_Std
		end # function NSE_θΨ

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NSE_WILMOT_θΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NSE_WILMOT_θΨ(hydro, N_Data, NiZ, optionₘ, θobs, Ψobs)

			NseWilmot_θΨ = zeros(Float64, NiZ)

			for iZ = 1:NiZ	
				θΨsim = zeros(Float64, N_Data[iZ])
				for iΨ = 1:N_Data[iZ]
					θΨsim[iΨ] = wrc.Ψ_2_θDual(optionₘ, Ψobs[iZ,iΨ], iZ, hydro)
				end
				NseWilmot_θΨ[iZ] = NSE_WILMOT( θobs[iZ,1:N_Data[iZ]], θΨsim[1:N_Data[iZ]])	
			end

			# Cumulating the objective function to get the overview
			NseWilmot_Mean = round(Statistics.mean(max.(NseWilmot_θΨ[1:NiZ] , 0.0)), digits=3)  # in case of negative value then it is set to 0
			NseWilmot_Std  = round(Statistics.std(max.(NseWilmot_θΨ[1:NiZ] , 0.0)), digits=3)   # in case of negative value then it is set to 0

		return NseWilmot_θΨ, NseWilmot_Mean, NseWilmot_Std
		end # function NSE_WILMOT_θΨ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NSE_KΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NSE_KΨ(hydro, N_Data, NiZ, optionₘ, Kobs, Ψobs)

			Nse_KΨ = zeros(Float64, NiZ)

			for iZ = 1:NiZ	
				KΨsim = zeros(Float64, N_Data[iZ])
				for iΨ = 1:N_Data[iZ]
					KΨsim[iΨ] = kunsat.Ψ_2_KUNSAT(optionₘ, Ψobs[iZ,iΨ], iZ, hydro)
				end
				Nse_KΨ[iZ] = 1.0 - NSE_MINIMIZE(log.(Kobs[iZ,1:N_Data[iZ]]), log.(KΨsim[1:N_Data[iZ]]))	
			end

			# Cumulating the objective function to get the overview
			Nse_Mean = round(Statistics.mean(max.(Nse_KΨ[1:NiZ] , 0.0)), digits=3)  # in case of negative value then it is set to 0
			Nse_Std  = round(Statistics.std(max.(Nse_KΨ[1:NiZ] , 0.0)), digits=3)   # in case of negative value then it is set to 0

		return Nse_KΨ, Nse_Mean, Nse_Std
		end # function NSE_KΨ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NSE_WILMOT_KΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NSE_WILMOT_KΨ(hydro, N_Data, NiZ, optionₘ, Kobs, Ψobs)

			NseWilmot_KΨ = zeros(Float64, NiZ)

			for iZ = 1:NiZ	
				KΨsim = zeros(Float64, N_Data[iZ])
				for iΨ = 1:N_Data[iZ]
					KΨsim[iΨ] = kunsat.Ψ_2_KUNSAT(optionₘ, Ψobs[iZ,iΨ], iZ, hydro)
				end
				NseWilmot_KΨ[iZ] = NSE_WILMOT(log.(Kobs[iZ,1:N_Data[iZ]]), log.(KΨsim[1:N_Data[iZ]]))	
			end

			# Cumulating the objective function to get the overview
			Nse_Mean = round(Statistics.mean(max.(NseWilmot_KΨ[1:NiZ] , 0.0)), digits=3)  # in case of negative value then it is set to 0
			Nse_Std  = round(Statistics.std(max.(NseWilmot_KΨ[1:NiZ] , 0.0)), digits=3)   # in case of negative value then it is set to 0

		return NseWilmot_KΨ, Nse_Mean, Nse_Std
		end # function NSE_WILMOT_KΨ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NSE_EFFICIENCY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NSE_EFFICIENCY(;Obs=Obs, Sim=Sim, Power=2.0)
			return Nse = 1 - NSE_MINIMIZE(Obs, Sim; Power=Power)
		end  # function: NSE_EFFICIENCY


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RELATIVE_ERR
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RELATIVE_ERR(; Obs=Obs, Sim=Sim)
			return Err = 1.0 - abs(Obs - Sim) / Obs
		end


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : LINEAR_REGRESSION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function LINEAR_REGRESSION(Xdata, Ydata)
			X = hcat(ones(length(Xdata)),Xdata)
			Y = Ydata
			return Intercept, Slope = inv(X'*X)*(X'*Y)
		end # function LINEAR_REGRESSION		function RELATIVE_ERR(;Obs=Obs, Sim=Sim)

end # module