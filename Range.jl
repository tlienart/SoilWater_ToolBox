module range

	# =============================================================
	#		MODULE: hydroparam
	# =============================================================
	module hydroparam
		θr_Min = 0.
		θr_Max = 0.25 # 0.2 or 0.25

		Ks_Mac_Max = 0.12 # 10.0 1.14[mm/s]

		# FEASIBLE RANGE OF KOSUGI HYDRAULIC PARAMETERS
			σMat_Min = 1.4 # 1.6
			σMat_Max = 4.3
		
			ΨkgMat_Min = 800.0 # [mm]
			ΨkgMat_Max = 60000.0 # 0.9 * 150000.0 #[mm]
		
			σMac_Min = 0.2
			σMac_Max = 0.8 #2.55
		
			ΨkgMac = 40. # 100 t0 10 [mm] 
			ΨkgMac_Min = 50. #[mm]
			ΨkgMac_Max = 390 #[mm]
	end  # module hydroparam
	# ............................................................

end