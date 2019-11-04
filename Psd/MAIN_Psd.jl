module psd
	import ..option, ..param, ..wrc, ..kunsat, ..cst, ..path, ..stats, ..psdFunc, ..psdInitiate, ..psdThetar
	using Statistics, BlackBoxOptim, JuliaDB

	# ======================================================================================
	#          PSD_START
	# ======================================================================================
	function PSD_MAIN(N_SoilSelect, Ψ_θΨ, θ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, Rpart, ∑Psd, N_Psd, hydro)

		# INITIATING THE PSD DATA		
			N_Psd, N_Psd_Max, Psd = psdInitiate.PSD_INITIATE(N_Psd, N_SoilSelect, ∑Psd)

			Nse_θr, θr_Psd = psd.psdThetar.MAIN_PSD_2_θr(N_SoilSelect, ∑Psd, hydro)

			println(θr_Psd)

				
	
			
		
				# ........................................................................






	end # function PSD_MAIN


end # module PSD