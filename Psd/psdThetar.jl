# =============================================================
#		MODULE: psdθr
# =============================================================
module psdThetar

	using ..param
	export PSD_2_θr

	# =========================================
	#       PSD -> θr 
	# =========================================
	# ∑Psd = 0.002 mm = Clay
		function PSD_2_θr(iSoil, ∑Psd; Psd_2_θr_α1 = param.psd.Psd_2_θr_α1, Psd_2_θr_α2 = param.psd.Psd_2_θr_α2)

			return θr_Psd = max( param.hydro.θr_Max * ( 1.0 - exp(- ( Psd_2_θr_α1 * (∑Psd[iSoil, param.psd.Psd_2_θr_Size] ^ Psd_2_θr_α2) ) ) ) , 0.0 )
		end # Function PSD_2_θr
	
end  # module name
# ............................................................