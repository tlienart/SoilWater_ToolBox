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

			Ψ_Max= 160000.0 #170000.0 # [mm] 160000.0 # min value is 150000 mm and oven dry would be the best value for the parameter 

		# Coeff_Φ_2_θs
			Coeff_Φ_2_θs = 0.98
			Coeff_θs_Max = 1.2

		# If constant
			ΨmMac = 40. # 100 t0 10 [mm]

			Ks_Min = 10.0 ^ -6.0  	# 0.000555556 wei [mm/s]
			Ks_Max = 2.0 # 10.0 ^ -4.0 	# 0.694444444 wei [mm/s]

			∇_θsMat_Min = 0.7

			σ_Min = 1.6 # 1.6
			σ_Max = 4.5

			Ψm_Min = 500.0 # [mm]
			Ψm_Max = 20000.0 # 0.9 * 150000.0 #[mm]

			∇_σMac_Min = 0.7 

			σMac_Min = 0.2
			σMac_Max = 0.8 #2.55
			∇_σMac_Max = 0.8

			ΨmMac_Min = 50. #[mm]
			ΨmMac_Max = 390 #[mm]

		# RELATIONSHIP BETWEEN PARAMETERS
			Pσ_1 = 0.5920
			Pσ_2 = 0.7679

	end  # module: hydro
	# ............................................................


	# =============================================================
	#		MODULE: infilt
	# =============================================================
	module infilt
		const SeIni_Output = [0.3 0.5 0.7] # [-] Different initial Se_Ini for plotting the infiltration curves 
		const Npoint_Infilt = 10 # Number of points for generating infiltration plots
		const Coeff_TransSteady = 5.0
	end  # module: infilt
	# ............................................................

end  # module param
# ............................................................