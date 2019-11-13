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
					Err_θr_Psd, θr_Psd, psdparam = OPTIMIZE_PSD_2_θr(N_SoilSelect, ∑Psd, hydro, psdparam)
		
				elseif option.psd.Psd_2_θr == "Cst" # <>=<>=<>=<>=<>
					θr_Psd = zeros(Float64, N_SoilSelect)
					for iSoil=1:N_SoilSelect
						psdparam.θr_Psd[iSoil] = param.psd.θr_Cst
					end
					fill!(psdparam.Psd_2_θr_α1, param.psd.θr_Cst) 
					fill!(psdparam.Psd_2_θr_α2, param.psd.θr_Cst)
					
				elseif option.psd.Psd_2_θr == "Param" # <>=<>=<>=<>=<>
					θr_Psd = zeros(Float64, N_SoilSelect)
					for iSoil=1:N_SoilSelect
						psdparam.θr_Psd[iSoil] = PSD_2_θr_FUNC(iSoil, ∑Psd)
					end
					# Putting the values Psd_2_θr_α1 & Psd_2_θr_α2 into psdparam
					fill!(psdparam.Psd_2_θr_α1, param.psd.Psd_2_θr_α1) 
					fill!(psdparam.Psd_2_θr_α2, param.psd.Psd_2_θr_α2)
				else # <>=<>=<>=<>=<>
					error("option.psd.Psd_2_θr = $option.psd.Psd_2_θr  not allowed option.psd.Psd_2_θr must be either (1)'Opt' or (2) 'Cst' or (3) 'Param' ")
				end # if option.psd.Psd_2_θr	

				return Err_θr_Psd, psdparam
			end # function PSD_2_θr(N_SoilSelect, ∑Psd, hydro)

		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>


		# =========================================
		#       PSD -> θr 
		# =========================================
		# ∑Psd = 0.002 mm = Clay
			function PSD_2_θr_FUNC(iSoil, ∑Psd; Psd_2_θr_α1 = param.psd.Psd_2_θr_α1, Psd_2_θr_α2 = param.psd.Psd_2_θr_α2)
				return θr_Psd = max( param.hydro.θr_Max * ( 1.0 - exp(- ( Psd_2_θr_α1 * (∑Psd[iSoil, param.psd.Psd_2_θr_Size] ^ Psd_2_θr_α2) ) ) ) , 0.0 )
			end # Function PSD_2_θr_FUNC

		
		# =========================================
		#       OPTIMIZE_PSD_2_θr 
		# =========================================
			function OPTIMIZE_PSD_2_θr(N_SoilSelect, ∑Psd, hydro, psdparam; Power=2)
				
				function OF(Psd_2_θr_α1, Psd_2_θr_α2, N_SoilSelect, ∑Psd, hydro)
					∑Rmse = 0.0
					for iSoil=1:N_SoilSelect
						θr_Psd = PSD_2_θr_FUNC(iSoil, ∑Psd; Psd_2_θr_α1=Psd_2_θr_α1, Psd_2_θr_α2=Psd_2_θr_α2)
						∑Rmse = abs(hydro.θr[iSoil] - θr_Psd) ^ Power
					end
				
					return ∑Rmse
				end # function OF =====

				Optimization = BlackBoxOptim.bboptimize(Param ->  OF(Param[1], Param[2], N_SoilSelect, ∑Psd, hydro) ; SearchRange =[ (param.psd.Psd_2_θr_α1_Min, param.psd.Psd_2_θr_α1_Max), (param.psd.Psd_2_θr_α2_Min, param.psd.Psd_2_θr_α2_Max)], NumDimensions=2, TraceMode=:silent)

				Psd_2_θr_α1 = BlackBoxOptim.best_candidate(Optimization)[1]
				Psd_2_θr_α2 = BlackBoxOptim.best_candidate(Optimization)[2]

				# Writing the values into psdparam
				fill!(psdparam.Psd_2_θr_α1, Psd_2_θr_α1) 
				fill!(psdparam.Psd_2_θr_α2, Psd_2_θr_α2)

				# COMPUTING THE OPTIMAL VALUE
					θr_Psd = zeros(Float64, N_SoilSelect)
					for iSoil=1:N_SoilSelect
						psdparam.θr_Psd[iSoil] = PSD_2_θr_FUNC(iSoil, ∑Psd; Psd_2_θr_α1=Psd_2_θr_α1, Psd_2_θr_α2=Psd_2_θr_α2)
					end

				# STATISTICS
					Nse_θr_Psd = stats.NASH_SUTCLIFFE_EFFICIENCY(;Obs=hydro.θr[1:N_SoilSelect], Sim=θr_Psd[1:N_SoilSelect])

					Err_θr_Psd = Array{Float64}(undef, (N_SoilSelect))
					for iSoil=1:N_SoilSelect
						Err_θr_Psd[iSoil] = stats.RELATIVE_ERR(;Obs=hydro.θr[iSoil], Sim=θr_Psd[iSoil])
					end
					
					println("    ~ Psd_2_θr_α1 = $(round(Psd_2_θr_α1,digits=3))  Psd_2_θr_α2 = $(round(Psd_2_θr_α2,digits=3)) ~")
					println("    ~ Nse_θr_Psd = $(round(Nse_θr_Psd,digits=3)) ~")

				return Err_θr_Psd, θr_Psd, psdparam
			end # function OPTIMIZE_PSD_2_θr
	

end  # module psdThetar
# ............................................................