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
		const Plot = false # <true>* plot or <false> no plotting
			const Plot_WaterRetentionCurve 	= true
			const Plot_BestLab 				= true
	
	# HYDRAULIC MODEL
		const HydroModel 		= "Kosugi" 		# <"Kosugi">* OR  <"Vangenuchten">
		const UnimodalBimodal 	= "Unimodal" 	# <"Unimodal" OR <"Bimodal>>

	# MODELS RUN
		const Id_Select 	= true 	# <true>* Select Id from the data OR <false> use all the data
		const θΨ 			= "Opt" # <"Opt">* Optimize hydraulic parameters from θ(ψ) OR <"File"> from save file OR <"No"> not available
		const KunsatΨ		= true 	#  <true>* Optimize hydraulic parameters from θ(ψ) & K(Ψ) OR <false>  
		const Psd 			= true 	# <true>* Derive θ(ψ) and/OR hydraulic parameters from Psd OR <false>
		const Infiltration 	= false 	# <true>* Derive θ(ψ) and/OR hydraulic parameters from Infiltration OR <false>
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
				const Model 		= "Chang2019Model" # <IMP> Intergranular Mixing Model OR <Chang2019Model> 
				const OptimizePsd 	= "Run" # <Single> or <All> or <Run>. <Single> =  optimize ξ1 & ξ2 for individual soils. <All> = derive universal parameters for all soils. <Run> = use parameters from Param.jl 
				const Psd_2_θr 		= "Param" # <Opt> optimises parameters α1 and α1; <Cst> uses θr = param.θr_Cst; <Param> uses α1 and α1 from parameters in Param.jl  # for new table model 1
				
				# For OptimizePsd = "Single"
					const ∑Psd_2_ξ1 = true  # <true> Use relationship between ξ1 and ∑Psd and do not optimize ξ1
					const ∑Psd_2_ξ2 = false # <true> Use relationship between ξ2 and ∑Psd and do not optimize ξ2
				
				# PLOTTING
				const Plotting = true # <true> Do plots as indicated below
					const Plot_σ_Ψkg = true
					const Plot_Kθ = true
					const Plot_IntergranularMixing = true
					const Plot_Pdf = true
					const Plot_∑Pdf = true
					const Plot_OFsingle_OFmeas = true
					const Plot_∑Psd_2_ξ2 = true
					const Plot_Ψkg_σmod = true
					const Plot_θr = true
							
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

end  # module option
# ............................................................