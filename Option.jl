"""  
	Recommended *options* are in *
"""
# =============================================================
#		MODULE: option
# =============================================================
module option
	# DOWNLAOD PACKAGES
		const DownloadPackage = false # <true> For first time user download packages required to run program OR <false>*

	# PLOTTING
		const Plot = true # <true>* plot or <false> no plotting
	
	# HYDRAULIC MODEL
		const HydroModel = 		"Kosugi" 		# <"Kosugi">* OR  <"Vangenuchten">
		const UnimodalBimodal = "Unimodal" 	# <"Unimodal" OR <"Bimodal>>

	# MODELS RUN
		const Id_Select 	= true 	# <true>* Select Id from the data OR <false> use all the data
		const θΨ 			= "Opt" # <"Opt">* Optimize hydraulic parameters from θ(ψ) OR <"File"> from save file OR <"No"> not available
		const KunsatΨ		= true 	#  <true>* Optimize hydraulic parameters from θ(ψ) & K(Ψ) OR <false>  
		const Psd 			= true 	# <true>* Derive θ(ψ) and/OR hydraulic parameters from Psd OR <false>
		const Infiltration 	= true 	# <true>* Derive θ(ψ) and/OR hydraulic parameters from Infiltration OR <false>
		const KsatModel 	= false # <true> Derive Ksat from θ(Ψ) OR <false>


		# =============================================================
		#		MODULE: hydroparam
		# =============================================================
			module hydro
			using ..option

				θsOpt 	= "Data" #  <Opt> Optimize θs OR <Data>* derived from Max θ(ψ) OR <Φ> which requires some correction
				if θsOpt == "Opt" && option.UnimodalBimodal == "Bimodal"
					println("\n NOT POSSIBLE: option.θsOpt == Opt && option.UnimodalBimodal = Bimodal")
					println("AUTO CORRECT: option.hydro.θsOpt = Data \n")
					θsOpt = "Data"
				end 

				θrOpt 	= "Psd" #  <"Opt">* Optimize θr OR  <"Psd"> Derived from particle size distribution: OR  θr=Cst <"Cst"> Fixed with value derived from param.hydro.θr

				KsOpt 	= "Data" #  <Opt> Optimize Ks (require KunsatΨ=true) OR <"Data"> from Max K(ψ)		
			end  # module hydroparam
		# ............................................................



		# =============================================================
		#		MODULE: psd
		# =============================================================
			module psd
				const Optimize 		= false # <true> Optimize the PSD model (require θΨ) OR <false>* derived from preset values 
				const HydroParam 	= false # <true> Optimize the hydraulic parameters from θ(ψ)psd OR <false>	
			end  # module psd
		# ............................................................



		# =============================================================
		#		MODULE: infiltration
		# =============================================================
			module inf
				const Optimize 		= false 		# <true> Optimize the parameters of the infiltration model (require θΨ) OR <false> derived from preset values
				const Model 		= "QuasiExact" 	# <"QuasiExact"> OR <"Simplified">*
				const HydroParam 	= false	 		# <true> Optimize the hydraulic parameters from θ(ψ)inf OR <false>
				const θini 			= false			# <true> derive the infiltration curves for different starting θini OR <false>
			end  # module inf
		# ............................................................
	


		# =============================================================
		#		MODULE: ksat
		# =============================================================
			module ksat
				const Optimize = true # <Opt> Optimize the parameters of Ks model (require KΨ) OR <false>* derived from preset values 
				const Ksat = false # <true> Derive Ksat from θ(Ψ) OR <false>
			end  # module ksat
			# ............................................................


	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	#	CHECKING OPTIONS
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	
				# KsOpt 
end  # module option
# ............................................................