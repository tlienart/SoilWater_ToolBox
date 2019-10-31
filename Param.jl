# =============================================================
#		MODULE: param
# =============================================================
module param
	N_iSoil_Simulations = 10 # maximum number of soils to be simulated (good for testing)

	# =============================================================
	#		MODULE:  θr FROM PSD
	# =============================================================
	module psd
		import ...option
		# WEIGHTED OF DERIVING THE LAST PARTICLE
			Wsubclay_Min = 0.1
			Wsubclay_Max = 1.0
			### optimized value: OptimizePsd = "Run" ###
			Subclay = 0.6934995359806453 

		# OPTIMISATION OF PSD
			Ψ_Max= 160000.0 #170000.0 # [mm] 160000.0 # min value is 150000 mm and oven dry would be the best value for the parameter 
			λ = 2.0 # exponent of the normalised Young-Laplace capillary equation # λ = 1 for new table model 1 ######

			ξ_Max = 3.0 

			ξ1_Min = 0.01 # 0.0 

			### optimized value: OptimizePsd = "Run" ###
				ξ1 = 9.040974907360946

			if option.psd.Model == "IMP"
				ξ1_Max =  20.0 
			elseif option.psd.Model == "Chang2019Model"
				ξ1_Max =  1.0 
			end # option.Chang2019Model
			
			ξ2_Min = 0.001   # 0.0 for new table model 1, 0.001 for new table model 4 ######
			ξ2_Max = 0.2     # 3.0 for new table model 1, 0.2   for new table model 4 ######

		# PEDOTRANSFERT FUNCTIONS
			# Relationship which computes θr from PSD
				Psd_2_θr_Size = 1  # size of particle size corresponding to clay fraction
				Psd_2_θr_α1_Min = 0.01
				Psd_2_θr_α1_Max = 100.0
				Psd_2_θr_α2_Min = 0.001
				Psd_2_θr_α2_Max = 10.0
				### optimized values: Psd_2_θr = "Param" ###
				Psd_2_θr_α1 = 16.01602133125399 # α1
				Psd_2_θr_α2 = 2.013125380534685 # α2 
		
			# If option.Psd_2_θr = false Residual θ kept constant during the model
				θr_Cst = 0.0 # 0.14515925    # is kept constant when not optimized
		
			# ξ1 constant
				P_ξ1 = 8.85 #3.514424626509076

		# Relationship which computes ξ2 from ∑Psd
			∑Psd_2_ξ2_Size = 2  # cumulative particle size fraction corresponding to very fine silt

			∑Psd_2_ξ2_β1_Min = 0.001 # for new table model 4   # ξ2_Min 
			∑Psd_2_ξ2_β1_Max = 0.1   # for new table model 4   # 1.0   

			∑Psd_2_ξ2_β2_Min  = 0.1 
			∑Psd_2_ξ2_β2_Max  = 5.0

			### optimized values: OptimizePsd = "Run" ###
			∑Psd_2_ξ2_β1 = 0.0874077451694647     
			∑Psd_2_ξ2_β2 = 0.9463302042937239		
	end  # module psd
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

			Ψ_MatrixMacro = 390 # 490[mm] determine when matrix and macro domain starts

			Ks_Min = 10.0 ^ -6.0  	# 0.000555556 wei [mm/s]
			Ks_Max = 0.7 # 10.0 ^ -4.0 	# 0.694444444 wei [mm/s]

			∇_θsMat_Min = 0.7

			σ_Min = 1.8 # 1.6
			σ_Max = 4.5

			Ψm_Min = 800.0 # [mm]
			Ψm_Max = 15000.0 # 0.9 * 150000.0 #[mm]

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
		const SeIni_Output = [0.0 0.2 0.4 0.6 0.8] # [-] Different initial Se_Ini for plotting the infiltration curves 
		const Npoint_Infilt = 300 # Number of points for generating infiltration plots
		const Coeff_TransSteady = 2.0
	end  # module: infilt
	# ............................................................

end  # module param
# ............................................................