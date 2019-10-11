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
		# θr = Cst
			θr = 0.0

		# Feasible range of Kosugi
			θr_Max = 0.25 # 0.2 or 0.25

		# Coeff_Φ_2_θs
			Coeff_Φ_2_θs = 0.98
			Coeff_θs_Max = 1.2

		# If constant
			ΨmMac = 40. # 100 t0 10 [mm]

			Ks_Min = 0.0001
			Ks_Max = 0.12 # 10.0 1.14[mm/s]
			∇_θs_Min = 0.7

			σ_Min = 1.4 # 1.6
			σ_Max = 4.3

			Ψm_Min = 800.0 # [mm]
			Ψm_Max = 60000.0 # 0.9 * 150000.0 #[mm]

			∇_σMac_Min = 0.7 

			σMac_Min = 0.2
			σMac_Max = 0.8 #2.55
			∇_σMac_Max = 0.8

			ΨmMac_Min = 50. #[mm]
			ΨmMac_Max = 390 #[mm]

	end  # module: hydro
	# ............................................................


	# =============================================================
	#		MODULE: infilt
	# =============================================================
	module infilt
		const SeIni_Output = [0.3 0.5 0.7] # [-] Different initial Se_Ini for plotting the infiltration curves 
		const Npoint_Infilt = 5 # Number of points for generating infiltration plots
		const Coeff_TransSteady = 5.0
	end  # module: infilt
	# ............................................................

end  # module param
# ............................................................