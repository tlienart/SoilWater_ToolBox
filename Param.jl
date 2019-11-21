# =============================================================
#		MODULE: param
# =============================================================
module param
	N_iSoil_Simulations = 10000 # maximum number of soils to be simulated (good for testing)

	# =============================================================
	#		MODULE:  θr FROM PSD
	# =============================================================
	module psd
		import ...option

			# OPTIMIZE θr: Relationship which computes θr from PSD
			θr_Cst = 0.0 # 0.14515925     # If option.Psd_2_θr = <false> then θr keep constant during the model

			Psd_2_θr_α1    = 16.01602133125399 # α1
				Psd_2_θr_α1_Min = 10.0
				Psd_2_θr_α1_Max = 100.0
			Psd_2_θr_α2    = 2.013125380534685 # α2  
				Psd_2_θr_α2_Min = 0.001
				Psd_2_θr_α2_Max = 10.0#5.0 we were using 10.0 in the IMP version!!!!!!!!!!!!

			Psd_2_θr_Size   = 1  # size of particle size corresponding to clay fraction

		# =============================================================
		#		MODULE: imp
		#		
		# =============================================================
		module imp
			# OPTIMISATION OF PSD
				Ψ_Max = 160000.0 # [mm] min value is 150000 mm and oven dry would be the best value  
				λ 	  = 2.0 	 # exponent of the normalised Young-Laplace capillary equation # λ = 1 for new table model 1 ######
				ξ_Max = 3.0 	 # ξ maximum physical value 

			# INTERGRANULAR MIXING PARTICLE MODEL
				ξ1 = 9.040974907360946
				ξ1_Min = 0.01 # 0.0 
				ξ1_Max =  20.0

				ξ2_Max = 0.2  

				∑Psd_2_ξ2_β1   = 0.0874077451694647 # Relationship which computes ξ2 from ∑Psd
					∑Psd_2_ξ2_β1_Min = 0.001 # for new table model 4   # ξ2_Min 
					∑Psd_2_ξ2_β1_Max = 0.1   # for new table model 4   # 1.0 
				∑Psd_2_ξ2_β2   = 0.9463302042937239
					∑Psd_2_ξ2_β2_Min = 0.1 
					∑Psd_2_ξ2_β2_Max = 5.0
				∑Psd_2_ξ2_Size = 2  # TODO cumulative particle size fraction corresponding to very fine silt
				Subclay        = 0.6934995359806453 # weighted of deriving the smallest particle size
					Subclay_Min = 0.1
					Subclay_Max = 1.0
		end  # module psi
		# ............................................................


		# =============================================================
		#		MODULE: chang
		# =============================================================
		module chan
			ξ1 = 0.5
				ξ1_Min =  0.0 # TODO ξ1_Min chang
				ξ1_Max =  1.0 # TODO ξ1_Max chang
		end  # module chan
		# ............................................................

	end  # module psd
	# ............................................................


	# =============================================================
	#		MODULE: hydro
	# =============================================================
	module hydro
		# θr = Cst
		θr = 0.0
			θr_Min = 0.0     # 0.2 or 0.25
			θr_Max = 0.25     # 0.2 or 0.25

		# Ks_Min = 10.0 ^ -6.0  	# 0.000555556 wei [mm/s]
		Ks_Max = 0.12 #0.7 # 10.0 ^ -4.0 	# 0.694444444 wei [mm/s]

		# Coeff_Φ_2_θs
			Coeff_Φ_2_θs = 0.95
			Coeff_θs_Max = 1.2 # Maximum value of θs when optimized

		# If constant
			Ψ_MatrixMacro = 390.0 # 490 [mm] determine when matrix and macro domain 
			Ψ_Max  = 160000.0 # [mm] min value is 150000 mm and oven dry would be the best value  


		# =============================================================
		#		MODULE: KOSUGI
		# =============================================================
		module kg
			σ_Min = 1.4 # 1.6
			σ_Max = 4.3

			σMac_Min = 0.2
			σMac_Max = 0.8 # 2.55
			
			Ψm_Min = 800.0   # [mm]
			Ψm_Max = 60000.0 # [mm] # 0.9*150000.0 
			
			ΨmMac  = 40.0  # 100 to 10 [mm]
				ΨmMac_Min = 50.0  # [mm]
				ΨmMac_Max = 390.0 # [mm]

			∇_θsMat_Min = 0.7
			∇_θsMat_Max = 1.0

			# RELATIONSHIP BETWEEN PARAMETERS
				Pσ_1 = 0.5920
				Pσ_2 = 0.7679

			
		end  # module kg
		# ............................................................


			
		# =============================================================
		#		MODULE: VAN GENUCHTEN
		# =============================================================
		module vg
			N_Min = 1.0001
			N_Max = 3.0

			Ψvg_Min = 10.0  # mm
			Ψvg_Max = 500000.0  # mm

		end  # module vg
		# ............................................................

	end  # module: hydro
	# ............................................................


	# =============================================================
	#		MODULE: infilt
	# =============================================================
	module infilt
		const SeIni_Output = [0.0 0.2 0.4 0.6 0.8] # [-] Different initial Se_Ini for plotting the infiltration curves 
		const Npoint_Infilt = 300 # Number of points for generating infiltration plots
		const Coeff_TransSteady = 2.0
	end  # module: infilt
	# ............................................................

end  # module param
# ............................................................