# =============================================================
#		MODULE: param
# =============================================================
module param

	# =============================================================
	#		MODULE:  θr FROM PSD
	# =============================================================
	module psd
		# Optimized values: Psd_2_θr = "Param"
		Psd_2_θr_Size = 1  # size of particle size corresponding to clay fraction
		Psd_2_θr_α1 = 16.01602133125399 # α1
		Psd_2_θr_α2 = 2.013125380534685 # α2 

		# Feasible range
			Psd_2_θr_α1_Min = 0.01
			Psd_2_θr_α1_Max = 100.0

			Psd_2_θr_α2_Min = 0.001
			Psd_2_θr_α2_Max = 10.0		
	end  # module thetarPsd
	# ............................................................


	# =============================================================
	#		MODULE: hydro
	# =============================================================
	module hydro
		# Feasible range of Kosugi
			θr_Min = 0.
			θr_Max = 0.25 # 0.2 or 0.25

			σMat_Min = 1.4 # 1.6
			σMat_Max = 4.3

			ΨkgMat_Min = 800.0 # [mm]
			ΨkgMat_Max = 60000.0 # 0.9 * 150000.0 #[mm]

			σMac_Min = 0.2
			σMac_Max = 0.8 #2.55

			ΨkgMac = 40. # 100 t0 10 [mm]

			ΨkgMac_Min = 50. #[mm]
			ΨkgMac_Max = 390 #[mm]

			Ks_Mac_Max = 0.12 # 10.0 1.14[mm/s]
		
	end  # module: hydro
	# ............................................................

end  # module param
# ............................................................