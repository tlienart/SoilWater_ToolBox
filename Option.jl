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
			const Plot_WaterRetentionCurve 	= true
			const Plot_BestLab 				= true
	
	# HYDRAULIC MODEL
		const HydroModel 		= "Kosugi" 		# <"Kosugi">* OR  <"Vangenuchten">
		const UnimodalBimodal 	= "Unimodal" 	# <"Unimodal" OR <"Bimodal>>

	# MODELS RUN
		const Id_Select 	= true 	# <true>* Select Id from the data OR <false> use all the data
		const θΨ 			= "Opt" # <"Opt">* Optimize hydraulic parameters from θ(ψ) OR <"File"> from save file OR <"No"> not available
		const KunsatΨ		= true 	#  <true>* Optimize hydraulic parameters from θ(ψ) & K(Ψ) OR <false>  
		const Psd 			= false 	# <true>* Derive θ(ψ) and/OR hydraulic parameters from Psd OR <false>
		const Infiltration 	= true 	# <true>* Derive θ(ψ) and/OR hydraulic parameters from Infiltration OR <false>
		const KsatModel 	= false # <true> Derive Ksat from θ(Ψ) OR <false>


		# =============================================================
		#		MODULE: hydroparam
		# =============================================================
			module hydro
			import ..option

				θsOpt 	= "Φ" #  <Opt> Optimize θs OR <Data>* derived from Max θ(ψ) OR <Φ> which requires some correction
				if θsOpt == "Opt" && option.UnimodalBimodal == "Bimodal"
					println("\n NOT POSSIBLE: option.θsOpt == Opt && option.UnimodalBimodal = Bimodal")
					println("AUTO CORRECT: option.hydro.θsOpt = Data \n")
					θsOpt = "Data"
				end 

				θrOpt 	= "Opt" #  <"Opt">* Optimize θr OR  <"Psd"> Derived from particle size distribution: OR  θr=Cst <"Cst"> Fixed with value derived from param.hydro.θr

				KsOpt 	= "Opt" #  <Opt> Optimize Ks (require KunsatΨ=true) OR <"Data"> from Max K(ψ)		
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
			module infiltration
				const OptimizeRun  		= "Run" # <"Opt">* Optimize hydraulic parameters from infiltration data <"Run"> run the inftration curves from known hydraulic parameters <"RunOpt"> run and optimise for comparison purposes
				const Model 			= "Simplified" 	# <"QuasiExact"> OR <"Simplified">*
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


		# =============================================================
		#		MODULE: psd
		# =============================================================
			module psd
				model = "IMP" # <IMP> Intergranular Mixing Model <Chang2019Model> 

				OptimizePsd = "All" # "Single" or "All" or "Run". "Single" =  optimize ξ1 & ξ2 for individual soils. "All" = derive universal parameters for all soils. "Run" = use parameters from Param.jl 
				
				Psd_2_θr = "Opt" # (1) "Opt" optimises parameters α1 and α1; (2) "Cst" uses θr = param.θr_Cst; (3) "Param" uses α1 and α1 from parameters in Param.jl  # for new table model 1

				# For OptimizePsd = "SINGLE"
					∑Psd_2_ξ1 = true # If "TRUE" we use relationship between ξ1 and ∑Psd and we do not optimize ξ1.

					∑Psd_2_ξ2 = false # If "TRUE" we use relationship between ξ2 and ∑Psd and we do not optimize ξ2.
				
				Psd_2_HydrauParam = true # If "TRUE" we optimize the hydraulic parameters from PSD
				
				Chang2019Model = false # when true computes WRC from PSD using Chang et al., 2019
				
				# PLOTTING ========================
					Plot_σ_Ψkg = true

					Plotting = true # Its plots
						Plot_Kθ = true
						Plot_IntergranularMixing = true
						Plot_Pdf = true
						Plot_∑Pdf = true
						Plot_OFsingle_OFmeas = true
						Plot_∑Psd_2_ξ2 = true
						Plot_Ψkg_σmod = true
						Plot_θr = true
								
					if OptimizePsd == "Single" 
						SubclayOpt = false
					else
						SubclayOpt = true # Determine if optimize an additional fraction < 2 μm clay content or if fixed deriving from a constant value param.Subclay
					end		
				end # module psd
			# ............................................................

end  # module option
# ............................................................