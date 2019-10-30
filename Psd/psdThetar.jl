# =============================================================
#		MODULE: psdθr
# =============================================================
module psdThetar

	using ..param

	export PSD_2_θr, OPTIMIZE_PSD_2_θr

	# =========================================
	#       PSD -> θr 
	# =========================================
	# ∑Psd = 0.002 mm = Clay
		function PSD_2_θr(iSoil, ∑Psd; Psd_2_θr_α1 = param.psd.Psd_2_θr_α1, Psd_2_θr_α2 = param.psd.Psd_2_θr_α2)

			return θr_Psd = max( param.hydro.θr_Max * ( 1.0 - exp(- ( Psd_2_θr_α1 * (∑Psd[iSoil, param.psd.Psd_2_θr_Size] ^ Psd_2_θr_α2) ) ) ) , 0.0 )
		end # Function PSD_2_θr

		
	# =========================================
	#       OPTIMIZE_PSD_2_θr 
	# =========================================
		function OPTIMIZE_PSD_2_θr(Nsample, θr, ∑Psd)
			θr_Psd = zeros(Float64, Nsample)

			function OF(Psd_2_θr_α1, Psd_2_θr_α2, Nsample, θr, Psd)
				Rmse = 0.0
				@simd for iSoil=1:Nsample
					θr_Psd[iSoil] = PSD_2_θr(∑Psd[iSoil]; Psd_2_θr_α1=Psd_2_θr_α1, Psd_2_θr_α2=Psd_2_θr_α2)
					Rmse += (θr_Psd[iSoil] - θr[iSoil]) ^ 2.0
				end
			
				return Rmse, θr_Psd
			end # function OF

			Optimization = BlackBoxOptim.bboptimize(Param ->  OF( Param[1],Param[2], Nsample, θr[1:Nsample], ∑Psd[1:Nsample, param.Psd_2_θr_Size])[1] ; SearchRange =[ (param.Psd_2_θr_α1_Min, param.Psd_2_θr_α1_Max), (param.Psd_2_θr_α2_Min, param.Psd_2_θr_α2_Max)], NumDimensions=2, TraceMode=:silent)

			Psd_2_θr_α1 = BlackBoxOptim.best_candidate(Optimization)[1]
			Psd_2_θr_α2 = BlackBoxOptim.best_candidate(Optimization)[2]

			~, θr_Psd = OF(Psd_2_θr_α1, Psd_2_θr_α2, Nsample, θr[1:Nsample], ∑Psd[1:Nsample, param.Psd_2_θr_Size])

			println("Optimize θr = Psd_2_θr_α1 = $Psd_2_θr_α1 ; Psd_2_θr_α2 = $Psd_2_θr_α2")
			return θr_Psd
		end # function OPTIMIZE_PSD_2_θr
	
end  # module psdThetar
# ............................................................


		
	