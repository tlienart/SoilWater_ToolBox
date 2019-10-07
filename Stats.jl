module stats
	using ..wrc, ..param
	export NASH_SUTCLIFFE_ERRORmin, NASH_SUTCLIFFE_θΨ, RELATIVEerr
	using Statistics

	function NASH_SUTCLIFFE_ERRORmin(Obs, Sim; Power=2.0)
		N = length(Obs)
		Obs_Mean = Statistics.mean(Obs[1:N])
		Obs_Mean_Err =  Statistics.sum(abs.(Obs_Mean .- Obs[1:N]) .^ Power)

		if Obs_Mean_Err < 0.000000000000001   
			Obs_Mean_Err = 1.0
		end
		
		Err = sum(abs.(Sim[1:N] - Obs[1:N]) .^ Power)

		return Of = (Err / Obs_Mean_Err ) 
	end # function NASH_SUTCLIFFE_ERRORmin


	function NASH_SUTCLIFFE_θΨ(Nsample, Nrpart, Ψ_Rpart, θ_Rpart, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac)
		Nse_Psd = zeros(Float64, Nsample)
		@simd for iSoil = 1:Nsample	
			θΨ = zeros(Float64, Nrpart[iSoil])
			@simd for iRpart = 1:Nrpart[iSoil]
				θΨ[iRpart] = wrc.kg.Ψ_2_θdual(Ψ_Rpart[iSoil,iRpart], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil])
			end
			Nse_Psd[iSoil] = 1.0 - stat.NASH_SUTCLIFFE_ERRORmin(θΨ[1:Nrpart[iSoil]], θ_Rpart[iSoil,1:Nrpart[iSoil]])	
		end
		# Cumulating the objective function to get the overview
		Nse_Psd_Mean = Statistics.mean(max.(Nse_Psd[1:Nsample],0.0))  # in case of negative value then it is set to 0
		Nse_Psd_Std  = Statistics.std(max.(Nse_Psd[1:Nsample],0.0))   # in case of negative value then it is set to 0
		#println( "Nse_Psd_Mean = $Nse_Psd_Mean, \n")
		println( "Nse_Psd_Mean_" * param.Name * " = $Nse_Psd_Mean,")
		println( "Nse_Psd_Std_" * param.Name * " = $Nse_Psd_Std, \n")
		return Nse_Psd, Nse_Psd_Mean
	end # function NASH_SUTCLIFFE_θΨ



	function RELATIVEerr(Obs, Sim)
		return Err = 1. - Statistics.abs(Obs - Sim ) / Obs
	end

end # module