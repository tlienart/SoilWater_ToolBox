# =============================================================
#		MODULE: option
# =============================================================
module option

        const Id_Select = true 	# <true>* Select Id from the data; <false> use all the data

		const θΨ        = "Opt" # <"Opt">* Optimize hydraulic parameters from θ(ψ); <"File"> from save file; <"No"> not available
      
        const Psd       = false	# <true>* Derive θ(ψ) AND/OR hydraulic parameters from Psd; <false>
		
        const Infilt    = true # <true>* Derive θ(ψ) AND/OR hydraulic parameters from Infiltration; <false>
     	
	# DOWNLAOD PACKAGES
        const DownloadPackage = false # <true> For first time user download packages required to run program; <false>*

	# PLOTTING
		const Plot = true # <true>* plot; <false> no plotting
			
	# =============================================================
	#		MODULE: hydro
	# =============================================================
		module hydro
			import ..option
		
			# HYDRAULIC MODEL
				const HydroModel      = "Kosugi" # <"Kosugi">*;  <"Vangenuchten">

				const UnimodalBimodal = "Unimodal" # <"Unimodal">; <"Bimodal>

				const θsOpt           = "Φ" #  <Opt> Optimize θs; <Φ> derived from total porosity which requires some correction from param.hydro.Coeff_Φ_2_θs; 
				
				const θrOpt           = "Opt" # <Opt> optimises; <Cst> uses θr = param.θr_Cst; <Param> uses α1 and α1 from parameters in Param.jl
				
				const KunsatΨ         = true #  <true>* Optimize hydraulic parameters from θ(ψ) & K(Ψ); <false>
					const KsOpt           = "Data" # <Opt> Optimize Ks (require KunsatΨ=true); <"Data"> derived from Max K(ψ)

			# PLOTTING
				const Plot_θΨ = true
			
			# if θsOpt == "Opt" && option.hydro.UnimodalBimodal == "Bimodal"
			# 	println("\n NOT POSSIBLE: option.θsOpt == Opt && option.hydro.UnimodalBimodal = Bimodal")
			# 	println("AUTO CORRECT: option.hydro.θsOpt = Data \n")
			# 	θsOpt = "Data"
			# end 
		end  # module hydro
		# ............................................................


		# =============================================================
		#		MODULE: psd      
		# =============================================================
		module psd
			const Model       = "IMP" # <IMP> Intergranular Mixing Model OR; <Chang2019Model>

			const OptimizePsd = "OptAllSoil" # <OptSingleSoil>; or <OptAllSoil>; or <Run>
			
			const Psd_2_θr    = "Opt" # <Opt> optimises parameters α1 and α1; <Cst> uses θr = param.θr_Cst; <Param> uses α1 and α1 from parameters in Param.jl 
			
			# FITTING THE PSD FUNCTION TO A HYDRAULIC MODEL			
				const HydroParam  = true # <true> Optimize the hydraulic parameters from θ(ψ)psd OR <false>
					const HydroModel      = "Kosugi" 		# <"Kosugi">*; OR  <"Vangenuchten">

					const UnimodalBimodal = "Unimodal" 	# <"Unimodal">; OR <"Bimodal>

                    const KunsatΨ = false 	# <false>  Not to change
                    const θsOpt   = "Known" # <Known>  Not to change
                    const θrOpt   = "Known" # <Known>  Not to change
                    const KsOpt   = "Known" # <Known>  Not to change

				# FOR OPTIMIZEPSD = "Single"
					const ∑Psd_2_ξ1 = true  # optimize ξ1
			
				# PLOTTING
					const Plot_Psd_θ_Ψ 	 = true # <true> include θ_Ψ values derived from IMP model or; <false> only θ_Ψ experimental values and fitted curve 
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
			# MODEL USED
				const Model            = "QuasiExact" 	# <"QuasiExact">; <"Best_Univ">

				const SingleDoubleRing = "Single"	# <"Double"> infiltration measured by single ring infiltrometer; <"Single"> infiltration measured by single ring infiltrometer
		
				const SortivityVersion = "NonInfinity" # <NonInfinity>; <Traditional>
				
				const SorptivityModel  = "Parlange" # <"Parlange"> strong non-linear diffusivity;  <"Crank"> constant diffusivity; <"Philip&Knight"> dirac delta-function diffusivity; <"Brutsaert"> moderate non-linear diffusivity,
				
				const Relation_σ_2_ψm  = false # <true> constrain param with Kosugi relation between σ and ψm or; <false>

				const Psd              = false # <true> constrain the opt hydraulic param with PSD; <false> ...

			# HYDRAULIC MODEL
				const HydroModel      = "Kosugi" 		# <"Kosugi">*; OR  <"Vangenuchten">
				const UnimodalBimodal = "Unimodal" 	# <"Unimodal">; OR <"Bimodal>
				const KunsatΨ         = false 	# <false>  Not to change
				const θsOpt           = "Known" # <Known>  Not to change
				const θrOpt           = "Known" # <Known>  Not to change
				const KsOpt           = "Known" # <Known>  Not to change
				
			# OUTPUT
				const OutputDimension = "3D" # <"1D"> such as by using single ring of <"3D"> by using double ring

				const OptimizeRun     = "Opt" # <"Opt">* Optimize hydraulic parameters from infiltration data <"Run"> run the inftration curves from known hydraulic parameters <"RunOptKs> run the inftration curves from known hydraulic parameters but optimize Ks only <"RunOpt"> run and optimise for comparison purposes <"RunOptKs"> run and optimise for comparison purposes without comparing Ks
			
			# PLOTTING
            const Plot_Sorptivity  = true # <true> or <false>
				
				const Plot_SeIni_Range = true # <true> computes infiltration curves for different SeIn set in param.infilt.SeIni_Output <false> no outputs
				
				const Plot_∑Infilt = true # <true> plots cumulative infiltration curves for experimental and derived data <false> no plots

				# =============================================================
				#		MODULE: besUniv
				# ==========================================================
		
				module bestUniv
					Continous = true # <true> this is the continous form derived by Joseph <false> Traditional non continous
				end  # module: besUniv

				# ............................................................

		end  # module infilt
		# ............................................................
	


		# =============================================================
		#		MODULE: ksat
		# =============================================================
			module ksat
				const Optimize = true # <Opt> Optimize the parameters of Ks model (require KΨ) OR <false>* derived from preset values
				
                const Ksat     = false # <true> Derive Ksat from θ(Ψ) OR <false>
			end  # module ksat
			# ............................................................

end  # module option
# ............................................................