# =============================================================
#		MODULE: option
# =============================================================
''' The recommended options are in * '''
module option

	# DOWNLAOD PACKAGES
		const DownloadPackage = false

	# HYDRAULIC MODEL
		const HydroModel = "vangenuchten" # [kosugi_dual]* or [kosugi_single]or [vangenuchten]

	# DATA AVAILABLE
		const θΨ = 		true  		# [true]* or [false]
		const KunsatΨ = false 		# [true]* or [false]
		const Psd = 	true 		# [true]* or [false]
		const InfT = 	true 		# [true]* or [false]
		const IdTrue = 	true 		# [true]* or [false]
		const OtherId = true 		# [true]* or [false]


		# =============================================================
		#		MODULE: hydroparam
		# =============================================================
			module hydroparam
				const Psd_2_θΨ = true
				const Psd_2_HydroParam = false		
			end  # module hydroparam
		# ............................................................


		# =============================================================
		#		MODULE: psd
		# =============================================================
			module psd
				const Psd_2_θΨ = true
				const Psd_2_HydroParam = false		
			end  # module psd
		# ............................................................



		# =============================================================
		#		MODULE: inf
		# =============================================================
			module inf
				const Model = "QuasiExact" # [QuasiExact] or [Simplified]
				const θΨ = false
				const HydroParam = false
				const Inf_θini = false	
			end  # module inf
		# ............................................................
	


		# =============================================================
		#		MODULE: ksat
		# =============================================================
			module ksat
				const Ksat = false # Derive Ksat from θ(Ψ)
			end  # module ksat
			# ............................................................

end  # module option
# ............................................................

