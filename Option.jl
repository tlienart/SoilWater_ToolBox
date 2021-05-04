# =============================================================
#		MODULE: global option
# ===========================================================
module option

		const HydroTranslateModel = false # <true>; <false>

      const Hypix       = false # <true>; <false>

      const Smap        = true # <true> ; <false>

      const BulkDensity = false # <true> <false>

      const θΨ          = :Opt # <:Opt>* Optimize hydraulic parameters from θ(Ψ); <:File> from save file; <:No> not available

      const Psd         = false	# <true> Derive θ(Ψ) AND/OR hydraulic parameters from Psd; <false>

      const Infilt      = false # <true> Derive θ(Ψ) AND/OR hydraulic parameters from Infiltration data; <false>

      const Temporary   = false # <true>; <false>

		const Jules = false # <true>; <false>
     	
	# DOWNLAOD PACKAGES
      const DownloadPackage = false # <true> For first time user download packages required to run program; <false>*

	# PLOTTING
      const Plot      = false # <true>* plot; <false> no plotting
      const Plot_Show = false # <true>* plot shown in VScode; <false>

	# =============================================================
	#		MODULE: global option
	# ===========================================================
		module data
			const Psd = true # <true> Particle size data is available; <false> Particle size data is available
		end
		
	# =============================================================
	#		MODULE: rock fragment
	# =============================================================
		module rockFragment
			const RockFragment = true # <true> make correction for rock fragment; <false> no correction for rock fragment
				const RockInjected = true # <true> rocks are injected in to the fine soils; <false> rocks are included in the bulk BulkDensity_Infilt
				
				const RockWettable = false # <true> rocks are wettable; <false> 
		end

	# =============================================================
	#		module: smap
	# =============================================================
		module smap
			const CorrectStone = true # <true> or <false>

			const CorrectStoneWetability = true # <true> or <false>
			
			const UsePointKosugiBimodal = true # <true> or <false>

			# AddPointKosugiBimodal = true # <true> or <false>

			AddPointKosugiBimodal = !(UsePointKosugiBimodal)

			const CombineData = true

			const Plot_Kunsat = false  # <true> or <false>
		end  # module: smap
		# ............................................................


	# =============================================================
	#		MODULE: hydro
	# =============================================================
		module hydro
			import ..option

			# HYDRAULIC MODEL
            const HydroModel      = :Kosugi # <:Kosugi>*; <:Vangenuchten>; <:BrooksCorey>; <:ClappHornberger>; <:VangenuchtenJules>

            const θrOpt           = :Opt # <:Opt> optimises; <:ParamPsd> derive from PSD and uses α1 and α1 from parameters in Param.jl; <:σ_2_θr>	
     
				const σ_2_Ψm          = :Constrained # <:Constrained> Ψm physical feasible range is computed from σ <:UniqueRelationship> Ψm is computed from σ; <:No> optimisation of σ & Ψm with no constraints

			# Min & Max from data
				const θs_MinFromData = false # <false> feasible range from GUI, <true> feasible range derive from data

				const Ks_MinMaxFromData = false # <false> feasible range from GUI, <true> feasible range derive from data

				
			# HAVE WE Kunsat(ψ)DATA
            const KunsatΨ         = true #  <true>* Optimize hydraulic parameters from θ(Ψ) & K(Ψ); <false>
					const KsOpt = :Opt # <:Opt> Optimize Ks (require KunsatΨ=true); <:Data> derived from Max K(Ψ)
					const Kunsat_JustRun = false

			# PLOTTING
				const Plot_θΨ = true
				const Plot_σ_Ψm = false
			
		end  # module hydro
		# ............................................................


	# =============================================================
	#		MODULE: psd      
	# =============================================================
		module psd
			const Model       = :IMP # <:IMP>* Intergranular Mixing Model; <:Chang2019Model>

			const OptimizePsd = :Run # <:OptSingleSoil>; <:OptAllSoil>; or <:Run>
			
			const Psd_2_θr    = :ParamPsd # <:Opt> optimises parameters α1 and α1; <:ParamPsd> uses α1 and α1 from parameters in Param.jl 

			# FOR OPTIMIZEPSD = :Single
					const ∑Psd_2_ξ1 = true  # optimize ξ1
			
			# FITTING THE PSD FUNCTION TO A HYDRAULIC MODEL			
				const HydroParam  = true # <true> Optimize the hydraulic parameters from θ(Ψ)psd OR <false>
				
				const HydroModel      = :Kosugi 		# <:Kosugi>*; <:Vangenuchten>; <:BrooksCorey>; <:ClappHornberger>

				const θsOpt           = :Φ #  <:Opt> Optimize θs; <:Φ> derived from total porosity which requires some correction from param.hydro.Coeff_Φ_2_θs;
				
            const θrOpt           = :ParamPsd # <:Opt> optimises; <:ParamPsd> derive from PSD and uses α1 and α1 from parameters in Param.jl; <:σ_2_θr>
				
            const σ_2_Ψm          = :Constrained # <:Constrained> Ψm physical feasible range is computed from σ  <:No> optimisation of σ & Ψm with no constraints
				
				const KunsatΨ         = false #  <true>* Optimize hydraulic parameters from θ(Ψ) & K(Ψ); <false>
					const KsOpt = :Opt # <:Opt> Optimize Ks (require KunsatΨ=true); <:Data> derived from Max K(Ψ)
					const Kunsat_JustRun = false

			# PLOTTING
            const Plot_Psd_θΨ    = true # <true>  plot θΨ of PSD; <false>
            const Plot_θr        = true #  <true>  plot θr data and model from Psd ; <false>	
				const Plot_IMP_Model = true # <true> ; plot IMP model results for publication; <false>
				
			# TABLE
				const Table_Psd_θΨ_θ = true # <true> derive θ values at prescribed Ψ
							
			# if OptimizePsd == :Single 
			# 	const SubclayOpt = false # Determine if optimize an additional fraction < 2 μm clay content or if fixed deriving from a constant value param.Subclay
			# else
			# 	const SubclayOpt = true # Determine if optimize an additional fraction < 2 μm clay content or if fixed deriving from a constant value param.Subclay
			# end			
		end  # module psd
		# ............................................................



		# =============================================================
		#		MODULE: infiltration
		# =============================================================
		module infilt
			# MODEL
            const DataSingleDoubleRing = :Single	# <:Double> infiltration measured by double ring infiltrometer; <:Single> infiltration measured by single ring infiltrometer

            const OptimizeRun          = :Opt # <:Opt>* Optimise hydraulic parameters from infiltration data; <:Run> run inftration curves from known hydraulic parameters; <:RunOptKs>  run inftration curves from known hydraulic parameters but optimize Ks only

            const Model                = :QuasiExact 	# <:QuasiExact> physical approach; <:Best_Univ> statistical improved approach

            const SortivityVersion     = :NonInfinity # <:NonInfinity> improved method; <:Traditional> old method with problems of infinity
				
            const SorptivitySplitModel = :Split # <:Split>; <:Split_η>
				
            const SorptivityModel      = :Parlange # <:Parlange> strong non-linear diffusivity;  <:Crank> constant diffusivity; <:Philip_Knight> dirac delta-function diffusivity; <:Brutsaert> moderate non-linear diffusivity,
				
			# DERIVING HYDRAULIC PARAMETERS FROM INFILTRATION TESTS
				const HydroModel      = :Kosugi # <:Kosugi>*; <:Vangenuchten>; <:BrooksCorey>; <:ClappHornberger>

				const θsOpt           = :Φ #  <:Opt> Optimize θs; <:Φ> derived from total porosity which requires some correction from param.hydro.Coeff_Φ_2_θs;
				
				const θrOpt           = :Opt # <:Opt> optimises; <:ParamPsd> derive from PSD and uses α1 and α1 from parameters in Param.jl; <:σ_2_θr>		

				const σ_2_Ψm          = :Constrained # <:Constrained> Ψm physical feasible range is computed from σ <:UniqueRelationship> Ψm is computed from σ; <:No> optimisation of σ & Ψm with no constraints
				
				const KunsatΨ         = true #  <true>* Optimize hydraulic parameters from θ(Ψ) & K(Ψ); <false>
					const KsOpt = :Opt # <:Opt> Optimize Ks (require KunsatΨ=true); <:Data> derived from Max K(Ψ)
					const Kunsat_JustRun = false
	
			# PLOTTING
            const Plot_Sorptivity       = true # <true> or <false>	
            const Plot_SeIni_Range      = true # <true> computes infiltration curves for different SeIn set in param.infilt.SeIni_Output <false> no outputs
            const Plot_∑Infilt          = true # <true> plots cumulative infiltration curves for experimental and derived data <false> no plots
            const Plot_θΨ               = true # <true>; <false>
            const Plot_Sorptivity_SeIni = true # <true> computes sorptivity curves as a function of Se <false> no outputs

				# =============================================================
				#		MODULE: besUniv
				# ==========================================================
				module bestFunc
					Continous = true # <true> this is the continous form derived by Joseph <false> Traditional non continous
				end  # module: besUniv

				# ............................................................

		end  # module infilt
		# ............................................................
	
	# =============================================================
	#		MODULE: hypix
	# =============================================================
	module hypix
		# DATA
			ClimateDataTimestep = "Daily" # <Hourly>; <Daily>

		# MODULES USED
         const RainfallInterception = true
         const Evaporation          = true
         const RootWaterUptake      = true
			const RootWaterUptakeComp  = true
			
		# SINK TERM 
         const LookupTable_Lai           = true # <false> Lai=constant; <true> Lai varies per month
         const LookUpTable_CropCoeficient = true # <false> CropCoeficient=constant; <true> CropCoeficient varies per month

		# HYDRAULIC MODEL 
			const θΨKmodel = :Kosugi # <:vanGenuchten>; <:Kosugi>

		# RICHARDS EQUATION
			const BottomBoundary = :Free # not working <:Free>; <:Pressure>
			
			const ∂R∂Ψ_Numerical = false # perform the derivatives numerically <true>; <false>

			# Adaptive time step
            # const Ψ_Constrain_K1     = false # <true>; <false>
				
            const AdaptiveTimeStep   = :Δθ # <:ΔΨ>; <:Δθ>
		
            const NormMin            = :Norm		#<:Norm>; <:Min>
				
            const Flag_ReRun         = true # <true>; <false> Rerun after updating the ΔT
				
            const Qbottom_Correction = true # <true> correction for the mass balance of the last cell
				
            # const NoConverge_Ψbest   = false # not working <true>; <false>* compute Q(Ψbest) when no convergence
		
		# RAINFALL INTERCEPTION MODEL
			const Lai_2_SintMax = false # <true> derive Sint_Sat from LAI_2_SINTMAX; <false> derive from input file

		# STEP WISE OPTIMIZATION
			const σ_2_Ψm = :No  # <:Constrained> Ψm physical feasible range is computed from σ <:UniqueRelationship> Ψm is computed from σ; <:No> optimisation of σ & Ψm with no constraints

			const σ_2_θr = true # <true> derive θr from σ <false>

			const θs_Opt = :No #  <:θs_Opt> θs is derived by multiplying a parameter to Max(θobs) for all profiles; <No>


		# CALIBRATION DATA AVAILABLE
         const calibr  = true # <true>; <false>

			const θobs_Hourly = true # θ data can be very large so we reduce the data to hourly

			const Signature_Run = false

		# TABLE true
         const Table = false 
            const Table_Discretization  = true
            const Table_Q               = true
            const Table_RootWaterUptake = true
            const Table_TimeSeries      = true
            const Table_Ψ               = true
            const Table_θ               = true
            const Table_TimeSeriesDaily = true
            const Tabule_θΨ             = true
            const Table_Climate         = true
			
		# PLOT OUTPUTS
         const Plot_Vegetation   = false
         const Plot_θΨK          = false
         const Plot_Interception = false
         const Plot_Other        = true
         const Plot_Sorptivity   = false
         const Plot_Hypix        = true
            const Plot_Climate      = true
            const Plot_θ            = true
            const Plot_Ψ            = true
            const Plot_Flux         = true
            const Plot_WaterBalance = true
            const Plot_ΔT           = true

	end  # module hypix
	# ............................................................

end  # module option
# ............................................................