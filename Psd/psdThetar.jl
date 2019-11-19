# =============================================================
#		MODULE: psdθr
# =============================================================
module psdThetar
	import ..param, ..stats, ..option
	import BlackBoxOptim
	export PSD_2_θr_FUNC, OPTIMIZE_PSD_2_θr

		# =========================================
		#       MAIN PSD -> θr 
		# =========================================
			function PSD_2_θr(N_SoilSelect, ∑Psd, hydro, psdparam)

				Err_θr_Psd = zeros(Float64, N_SoilSelect)

				if option.psd.Psd_2_θr == "Opt" && option.θΨ ≠ "No"
					psdparam = OPTIMIZE_PSD_2_θr(N_SoilSelect, ∑Psd, hydro, psdparam)
		
				elseif option.psd.Psd_2_θr == "Cst" # <>=<>=<>=<>=<>
					θr_Psd =  Array{Float64}(undef, N_SoilSelect)
					fill!(psdparam.θr_Psd, param.psd.θr_Cst)
					fill!(psdparam.Psd_2_θr_α1, 0.0) 
					fill!(psdparam.Psd_2_θr_α2, 0.0)
					
				elseif option.psd.Psd_2_θr == "Param" # <>=<>=<>=<>=<>
					θr_Psd =  Array{Float64}(undef, N_SoilSelect)
					for iSoil=1:N_SoilSelect
						psdparam.θr_Psd[iSoil] = PSD_2_θr_FUNC(iSoil, ∑Psd)
					end
					# Putting the values Psd_2_θr_α1 & Psd_2_θr_α2 into psdparam
					fill!(psdparam.Psd_2_θr_α1, param.psd.Psd_2_θr_α1) 
					fill!(psdparam.Psd_2_θr_α2, param.psd.Psd_2_θr_α2)

					# STATISTICS
					if option.θΨ ≠ "No"
						Nse_θr_Psd = stats.NASH_SUTCLIFFE_EFFICIENCY(;Obs=hydro.θr[1:N_SoilSelect], Sim=psdparam.θr_Psd[1:N_SoilSelect])

						println("    ~ Nse_θr_Psd = $(round(Nse_θr_Psd,digits=3)) ~")

						for iSoil=1:N_SoilSelect
							psdparam.Err_θr_Psd[iSoil] = stats.RELATIVE_ERR(;Obs=hydro.θr[iSoil], Sim=psdparam.θr_Psd[iSoil])
						end
					end
				
				else # <>=<>=<>=<>=<>
					error("option.psd.Psd_2_θr = $option.psd.Psd_2_θr  not allowed option.psd.Psd_2_θr must be either (1)'Opt' or (2) 'Cst' or (3) 'Param' ")
				end # if option.psd.Psd_2_θr	

				return psdparam
			end # function PSD_2_θr(N_SoilSelect, ∑Psd, hydro)

		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>


		# =========================================
		#       PSD -> θr 
		# =========================================
			function PSD_2_θr_FUNC(iSoil, ∑Psd; Psd_2_θr_Size=param.psd.Psd_2_θr_Size, θr_Max=param.hydro.θr_Max, Psd_2_θr_α1=param.psd.Psd_2_θr_α1, Psd_2_θr_α2=param.psd.Psd_2_θr_α2)

				return θr_Psd = max(θr_Max * ( 1.0 - exp(- ( Psd_2_θr_α1 * (∑Psd[iSoil,Psd_2_θr_Size] ^ Psd_2_θr_α2) ) ) ) , 0.0 )
				
			end # Function PSD_2_θr_FUNC

		
		# =========================================
		#       OPTIMIZE_PSD_2_θr 
		# =========================================
			function OPTIMIZE_PSD_2_θr(N_SoilSelect, ∑Psd, hydro, psdparam; Power=4)
				
				function OF(Psd_2_θr_α1, Psd_2_θr_α2, N_SoilSelect, ∑Psd, hydro)
					∑Rmse = 0.0
					for iSoil=1:N_SoilSelect
						θr_Psd = PSD_2_θr_FUNC(iSoil, ∑Psd; Psd_2_θr_α1=Psd_2_θr_α1, Psd_2_θr_α2=Psd_2_θr_α2)

						∑Rmse = abs(hydro.θr[iSoil] - θr_Psd) ^ Power
					end
				
					return ∑Rmse
				end # function OF =====

				SearchRange = [ (param.psd.Psd_2_θr_α1_Min, param.psd.Psd_2_θr_α1_Max), (param.psd.Psd_2_θr_α2_Min, param.psd.Psd_2_θr_α2_Max)]

				Optimization = BlackBoxOptim.bboptimize(Param ->  OF(Param[1], Param[2], N_SoilSelect, ∑Psd, hydro) ; SearchRange=SearchRange, NumDimensions=2, TraceMode=:silent)

				Psd_2_θr_α1 = BlackBoxOptim.best_candidate(Optimization)[1]
				Psd_2_θr_α2 = BlackBoxOptim.best_candidate(Optimization)[2]

				# Writing the values into psdparam
				fill!(psdparam.Psd_2_θr_α1, Psd_2_θr_α1) 
				fill!(psdparam.Psd_2_θr_α2, Psd_2_θr_α2)

				# COMPUTING THE OPTIMAL VALUE
					for iSoil=1:N_SoilSelect
						psdparam.θr_Psd[iSoil] = PSD_2_θr_FUNC(iSoil, ∑Psd; Psd_2_θr_α1=psdparam.Psd_2_θr_α1[iSoil], Psd_2_θr_α2=psdparam.Psd_2_θr_α2[iSoil])
					end

				# STATISTICS
					Nse_θr_Psd = stats.NASH_SUTCLIFFE_EFFICIENCY(;Obs=hydro.θr[1:N_SoilSelect], Sim=psdparam.θr_Psd[1:N_SoilSelect])

					for iSoil=1:N_SoilSelect
						psdparam.Err_θr_Psd[iSoil] = stats.RELATIVE_ERR(;Obs=hydro.θr[iSoil], Sim=psdparam.θr_Psd[iSoil])
					end
					
					println("    ~ Psd_2_θr_α1 = $(round(Psd_2_θr_α1,digits=3)) ;  Psd_2_θr_α2 = $(round(Psd_2_θr_α2,digits=3)) ;  Nse_θr_Psd = $(round(Nse_θr_Psd,digits=3)) ~")

				return psdparam
			end # function OPTIMIZE_PSD_2_θr
	
end  # module psdThetar
# ............................................................