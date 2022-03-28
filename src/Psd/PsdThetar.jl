# =============================================================
#		MODULE: psdθr
# =============================================================
module psdThetar
	import ..stats
	import BlackBoxOptim
	export PSD_2_θr_FUNC, OPTIMIZE_PSD_2_θr, PSD_2_θr

		# =========================================
		#       MAIN PSD -> θr 
		# =========================================
			function PSD_2_θr(∑Psd, hydro, hydroPsd, NiZ, option, param, paramPsd)

				Err_θr_Psd = zeros(Float64, NiZ)

				if option.psd.Psd_2_θr⍰=="Opt" && option.run.HydroLabθΨ⍰ ≠ "No"
					paramPsd = OPTIMIZE_PSD_2_θr(∑Psd, hydro, hydroPsd, NiZ, param, paramPsd)
							
				elseif option.psd.Psd_2_θr⍰ == "ParamPsd" # <>=<>=<>=<>=<>
					θr_Psd =  fill(0.0::Float64, NiZ)
					for iZ=1:NiZ
						paramPsd.θr_Psd[iZ] = PSD_2_θr_FUNC(∑Psd, hydroPsd, iZ, param)
					end
					# Putting the values Psd_2_θr_α1 & Psd_2_θr_α2 into paramPsd
					fill!(paramPsd.Psd_2_θr_α1, param.psd.Psd_2_θr_α1) 
					fill!(paramPsd.Psd_2_θr_α2, param.psd.Psd_2_θr_α2)
				
				else # if θr is not being optimised <>=<>=<>=<>=<>				
					θr_Psd =  fill(0.0::Float64, NiZ)
					fill!(paramPsd.θr_Psd, param.hydroPsd.θr)
					fill!(paramPsd.Psd_2_θr_α1, 0.0) 
					fill!(paramPsd.Psd_2_θr_α2, 0.0)
					for iZ=1:NiZ
						 paramPsd.θr_Psd[iZ] = hydroPsd.θr_Psd[iZ] 
					end
				end # if option.psd.Psd_2_θr⍰
				
				# TODO: Needs to hormonize hydroPsd with paramPsd
					for iZ=1:NiZ
						hydroPsd.θr[iZ] = paramPsd.θr_Psd[iZ] 
					end

				# STATISTICS
					if option.psd.Psd_2_θr⍰=="Opt" && option.run.HydroLabθΨ⍰ ≠ "No"
						Nse_θr_Psd = stats.NSE_EFFICIENCY(;Obs=hydro.θr[1:NiZ], Sim=paramPsd.θr_Psd[1:NiZ])

						for iZ=1:NiZ
							paramPsd.Err_θr_Psd[iZ] = stats.RELATIVE_ERR(;Obs=hydro.θr[iZ], Sim=paramPsd.θr_Psd[iZ])
						end

						println("    	~ Nse_θr_Psd = $(round(Nse_θr_Psd,digits=3)) ~ \n")
					end
					
			return hydroPsd, paramPsd
			end # function PSD_2_θr(NiZ, ∑Psd, hydroPsd)

		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>


		# =========================================
		#       PSD -> θr 
		# =========================================
			function PSD_2_θr_FUNC(∑Psd, hydroPsd, iZ, param; Psd_2_θr_Size=param.psd.Psd_2_θr_Size, Psd_2_θr_α1=param.psd.Psd_2_θr_α1, Psd_2_θr_α2=param.psd.Psd_2_θr_α2)
				return θr_Psd = max(hydroPsd.θr_Max[iZ] * (1.0 - exp(- (Psd_2_θr_α1 * (∑Psd[iZ,Psd_2_θr_Size] ^ Psd_2_θr_α2) ) )), 0.0)		
			end # Function PSD_2_θr_FUNC

		
		# =========================================
		#       OPTIMIZE_PSD_2_θr 
		# =========================================
			function OPTIMIZE_PSD_2_θr(∑Psd, hydro, hydroPsd, NiZ, param, paramPsd; Power=2)
				
				function OF(Psd_2_θr_α1, Psd_2_θr_α2, ∑Psd, hydro, NiZ, param, paramPsd)
					∑Rmse = 0.0
					for iZ=1:NiZ
						paramPsd.θr_Psd[iZ] = PSD_2_θr_FUNC(∑Psd, hydroPsd, iZ, param; Psd_2_θr_α1=Psd_2_θr_α1, Psd_2_θr_α2=Psd_2_θr_α2)

						∑Rmse += abs(hydro.θr[iZ] - paramPsd.θr_Psd[iZ]) ^ Power
					end
				
				return ∑Rmse
				end # function OF ======================================

				SearchRange = [(param.psd.Psd_2_θr_α1_Min, param.psd.Psd_2_θr_α1_Max), (param.psd.Psd_2_θr_α2_Min, param.psd.Psd_2_θr_α2_Max)]

				Optimization = BlackBoxOptim.bboptimize(Param -> OF(Param[1], Param[2], ∑Psd, hydro, NiZ, param, paramPsd) ; SearchRange=SearchRange, NumDimensions=2, TraceMode=:silent)

				Psd_2_θr_α1 = BlackBoxOptim.best_candidate(Optimization)[1]
				Psd_2_θr_α2 = BlackBoxOptim.best_candidate(Optimization)[2]

				# Writing the values into paramPsd
				fill!(paramPsd.Psd_2_θr_α1, Psd_2_θr_α1) 
				fill!(paramPsd.Psd_2_θr_α2, Psd_2_θr_α2)

				# COMPUTING THE OPTIMAL VALUE
					for iZ=1:NiZ
						paramPsd.θr_Psd[iZ] = PSD_2_θr_FUNC(∑Psd, hydroPsd, iZ, param; Psd_2_θr_α1=paramPsd.Psd_2_θr_α1[iZ], Psd_2_θr_α2=paramPsd.Psd_2_θr_α2[iZ])
					end
					println("    == Optimizing θr from PSD ==")
					println("    	~ Psd_2_θr_α1 = $(round(Psd_2_θr_α1,digits=3)) ;  Psd_2_θr_α2 = $(round(Psd_2_θr_α2,digits=3)) ~")

			return paramPsd
			end # function OPTIMIZE_PSD_2_θr
	
end  # module psdThetar
# ............................................................