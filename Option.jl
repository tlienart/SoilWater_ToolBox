# =============================================================
#		MODULE: option
# =============================================================
module option
	# MODELS RUN
        const Id_Select = true 	# <true>* Select Id from the data OR <false> use all the data
        const θΨ        = "Opt" # <"Opt">* Optimize hydraulic parameters from θ(ψ) OR <"File"> from save file OR <"No"> not available
      
        const Psd       = true 	# <true>* Derive θ(ψ) and/OR hydraulic parameters from Psd OR <false>
        const Infilt    = false 	# <true>* Derive θ(ψ) and/OR hydraulic parameters from Infiltration OR <false>
     	
	# DOWNLAOD PACKAGES
		const DownloadPackage = false # <true> For first time user download packages required to run program OR <false>*

	# PLOTTING
		const Plot = true # <true>* plot or <false> no plotting

	
			
		# =============================================================
		#		MODULE: hydro
		# =============================================================
		module hydro
			import ..option
		
		# HYDRAULIC MODEL
            const HydroModel      = "Kosugi" 		# <"Kosugi">* OR  <"Vangenuchten">

            const UnimodalBimodal = "Unimodal" 	# <"Unimodal" OR <"Bimodal>

            const KunsatΨ         = true 	#  <true>* Optimize hydraulic parameters from θ(ψ) & K(Ψ) OR <false>

            const θsOpt           = "Opt" #  <Opt> Optimize θs OR <Data>* derived from Max θ(ψ) OR <Φ> which requires some correction from param.hydro.Coeff_Φ_2_θs

            const θrOpt           = "Opt"  # <Opt> optimises; <Cst> uses θr = param.θr_Cst; <Psd> uses α1 and α1 from parameters in Param.jl

            const KsOpt           = "Opt" #  <Opt> Optimize Ks (require KunsatΨ=true) OR <"Data"> from Max K(ψ)

			# PLOTTING
				const Plot_θΨ = true
			
			if θsOpt == "Opt" && option.hydro.UnimodalBimodal == "Bimodal"
				println("\n NOT POSSIBLE: option.θsOpt == Opt && option.hydro.UnimodalBimodal = Bimodal")
				println("AUTO CORRECT: option.hydro.θsOpt = Data \n")
				θsOpt = "Data"
			end 
		end  # module hydro
		# ............................................................



		# =============================================================
		#		MODULE: psd      
		# =============================================================
		module psd
			const Model       = "IMP" # <IMP> Intergranular Mixing Model OR <Chang2019Model>
			const OptimizePsd = "OptAllSoil" # <OptSingleSoil> or <OptAllSoil> or <Run>
			const Psd_2_θr    = "Opt" # <Opt> optimises parameters α1 and α1; <Cst> uses θr = param.θr_Cst; <Param> uses α1 and α1 from parameters in Param.jl 
			
			# FITTING THE PSD FUNCTION TO A HYDRAULIC MODEL
			const HydroParam  = false # <true> Optimize the hydraulic parameters from θ(ψ)psd OR <false>
				const HydroModel      = "Kosugi" 		# <"Kosugi">* OR  <"Vangenuchten">

				const UnimodalBimodal = "Unimodal" 	# <"Unimodal" OR <"Bimodal>

				const KunsatΨ         = false 	#  <true>* Optimize hydraulic parameters from θ(ψ) & K(Ψ) OR <false>

				const θsOpt           = "Φ" #  <Opt> Optimize θs OR <Data>* derived from Max θ(ψ) OR <Φ> which requires some correction from param.hydro.Coeff_Φ_2_θs

				const θrOpt           = "Psd"  # <Opt> optimises; <Cst> uses θr = param.θr_Cst; <Psd> uses α1 and α1 from parameters in Param.jl

				const KsOpt           = "Data" #  <Opt> Optimize Ks (require KunsatΨ=true) OR <"Data"> from Max K(ψ)

			# FOR OPTIMIZEPSD = "Single"
				const ∑Psd_2_ξ1 = true  # optimize ξ1
		
			# PLOTTING
				const Plot_Psd_θ_Ψ 	 = true # <true> include θ_Ψ values derived from IMP model or <false> only θ_Ψ experimental values and fitted curve 
				const Plot_θr 	   	 = true # plot θr data and model from Psd 
				const Plot_IMP_model = true # plot IMP model results for publication
						
			if OptimizePsd == "Single" 
				const SubclayOpt = false # Determine if optimize an additional fraction < 2 μm clay content or if fixed deriving from a constant value param.Subclay
			else
				const SubclayOpt = true # Determine if optimize an additional fraction < 2 μm clay content or if fixed deriving from a constant value param.Subclay
			end			
		end  # module psd
		# ............................................................



		# =============================================================
		#		MODULE: infiltration
		# =============================================================
		module infilt
			const OptimizeRun  		= "Run" # <"Opt">* Optimize hydraulic parameters from infiltration data <"Run"> run the inftration curves from known hydraulic parameters <"RunOpt"> run and optimise for comparison purposes
			const Model 			= "Best_JJ" 	# <"QuasiExact"> OR <"Best_JJ">*
			const HydroParam 		= false	 		# <true> Optimize the hydraulic parameters from θ(ψ)inf OR <false>
			const Dimension			= "3D"	# <"3D"> infiltration rate by using single ring infiltrometer <"1D"> making 1D infiltration rate by using double ring infiltrometer
			const Relation_σ_2_ψm   = false # <true> one parameter will be optimized only
			const SeIni_Range 		= true # <true> computes infiltration curves for different SeIn set in param.infilt.SeIni_Output <false> no outputs
			const Sorptivity		= "Parlange" # <Parlange> or <>
		end  # module infilt
		# ............................................................
	


		# =============================================================
		#		MODULE: ksat
		# =============================================================
			module ksat
				const Optimize = true # <Opt> Optimize the parameters of Ks model (require KΨ) OR <false>* derived from preset values 
				const Ksat = false # <true> Derive Ksat from θ(Ψ) OR <false>
			end  # module ksat
			# ............................................................

end  # module option
# ............................................................