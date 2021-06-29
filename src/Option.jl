# =============================================================
#		module option
# ===========================================================
module options
	# What available data we have?
	struct DATA
		HydroParam::Bool
		Infilt::Bool
		Kθ::Bool
		Pedological⍰::Symbol
		Psd::Bool
		RockWetability::Bool
		SimulationKosugiθΨK::Bool
		θΨ::Bool
		Φ⍰::Symbol
	end # struct DATA

	# What model wanting to run
	mutable struct RUN
      ChangeHydroModel::Bool
      HydroLabθΨ⍰::Symbol
      Hypix::Bool
      Infilt::Bool
      IntergranularMixingPsd ::Bool
      RockCorection::Bool
		Smap::Bool
		Smap2Hypix::Bool
      Temporary::Bool
	end

	# How to input the data
	mutable struct DATAFROM
		Jules::Bool
	end
	mutable struct OTHER
		DownloadPackage::Bool
		Ploting::Bool
		PlotVscode::Bool
		DataPrivateShare::String
	end
	mutable struct SMAP
		Plot_Kunsat
	end
	mutable struct HYDRO
		HydroModel⍰::Symbol
		HydroModel_List::Vector{Symbol}
		θsOpt⍰::Symbol
		θrOpt⍰::Symbol
		σ_2_Ψm⍰::Symbol
		Plot_θΨ::Bool
		Plot_σ_Ψm::Bool
	end
	mutable struct PSD
		Model⍰::Symbol
		OptimizePsd⍰::Symbol
		Psd_2_θr⍰::Symbol
		∑Psd_2_ξ1::Bool
		HydroModel⍰::Symbol
		θsOpt⍰::Symbol
		θrOpt⍰::Symbol
		σ_2_Ψm⍰::Symbol
		Plot_Psd_θΨ::Bool
		Plot_θr::Bool
		Plot_IMP_Model::Bool
		Table_Psd_θΨ_θ::Bool
	end
	mutable struct ROCKFRAGMENT
		CorectStoneRockWetability::Bool
		RockInjectedIncluded⍰::Symbol
	end
	mutable struct INFILT
		DataSingleDoubleRing⍰::Symbol  
		OptimizeRun⍰::Symbol  
		Model⍰::Symbol
		BestUniv_Continous::Bool                	  
		SorptivitySplitModel⍰::Symbol  
		SorptivityModel⍰::Symbol     
		HydroModel⍰::Symbol       
		θsOpt⍰::Symbol            
		θrOpt⍰::Symbol            
		σ_2_Ψm⍰::Symbol               
		Plot_Sorptivity::Bool        	
		Plot_SeIni_Range::Bool       
		Plot_∑Infilt::Bool           
		Plot_θΨ::Bool                
		Plot_Sorptivity_SeIni::Bool  
	end
	mutable struct HYPIX
		ClimateDataTimestep⍰
		RainfallInterception
		Evaporation
		RootWaterUptake
		RootWaterUptakeComp
		LookupTable_Lai
		LookUpTable_CropCoeficient
		HydroModel⍰
		BottomBoundary
		∂R∂Ψ_Numerical
		AdaptiveTimeStep
		NormMin
		Flag_ReRun
		Qbottom_Correction
		Lai_2_SintMax
		σ_2_Ψm⍰
		σ_2_θr
		θs_Opt
		Optimisation::Bool
		θobs::Bool
		θobs_Average::Bool
		θobs_Hourly::Bool
		Signature_Run
		Table::Bool
		Table_Discretization::Bool
		Table_Q::Bool
		Table_RootWaterUptake::Bool
		Table_TimeSeries::Bool
		Table_Ψ::Bool
		Table_θ::Bool
		Table_TimeSeriesDaily::Bool
		Tabule_θΨ::Bool
		Table_Climate::Bool
		Plot_Vegetation::Bool
		Plot_θΨK::Bool
		Plot_Interception::Bool
		Plot_Other::Bool
		Plot_Sorptivity::Bool
		Plot_Hypix::Bool
		Plot_Climate::Bool
		Plot_θ::Bool
		Plot_Ψ::Bool
		Plot_Flux::Bool
		Plot_WaterBalance::Bool
		Plot_ΔT::Bool
	end
		mutable struct OPTION
			data::DATA
			dataFrom::DATAFROM
			hydro::HYDRO
			hyPix::HYPIX
			infilt::INFILT
			other::OTHER
			psd::PSD
			rockFragment::ROCKFRAGMENT
			run::RUN
			smap::SMAP
		end # struct OPTION
	
	#__________________________________________________________________
	#..................................................................

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPTION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OPTIONS()
			# =============================================================
			# 		OTHER
			# =============================================================
				# Download packages
					DownloadPackage = false # <true> For first time user download packages required to run program; <false>*

				# Plotting
					Ploting   = true # <true>* plot; <false> no plotting
					PlotVscode = true # <true>* plot shown in VScode; <false>
					DataPrivateShare = "Private" # <"Private"> data kept private in GitHub; <"Share"> date share in GitHub

			other = OTHER(DownloadPackage, Ploting, PlotVscode, DataPrivateShare)

			# =============================================================
			# 		DATA
			#      What data do we have ?
			# =============================================================
            HydroParam          = false # <true> ; <false>
            Infilt              = true # <true> ; <false>
            Jules               = false # <true> or <false>
            Kθ                  = true # <true> ; <false>
            Pedological⍰        = :Core # <:Core>; from traditional data set <:Smap> from Smap; <No> no data available
            Psd                 = true # <true> ; <false>
            RockWetability      = true # <true> ; <false>
            SimulationKosugiθΨK = false # <true> or <false>
            Φ⍰                  = :ρᵦ # <:ρᵦ> BulkDensity data; <:Φ> TotalPorosity; <:No> no data
            θΨ                  = true # <true> ; <false>

			data = DATA(HydroParam, Infilt, Kθ , Pedological⍰, Psd, RockWetability, SimulationKosugiθΨK, θΨ,  Φ⍰)


			# =============================================================
			# 		DATA FROM
			#      How to read data ?
			# =============================================================
				Jules = false

			dataFrom = DATAFROM(Jules)

			# =============================================================
			# 		RUN
			#      What model wanting to run ?
			# =============================================================
            ChangeHydroModel       = false # Not yet operational<true>; <false>
            HydroLabθΨ⍰           = :Opt # <:Opt>* Optimize hydraulic parameters from θ(Ψ); <:File> from save file; <:Run> just run <:No> not available
            Hypix                  = false # <true>; <false>
            Infilt                 = false # <true>; <false>
            IntergranularMixingPsd = true # <true>; <false>
            RockCorection          = false # <true> make correction for rock fragment; <false> no correction for rock fragment
            Smap                   = false  # <true>; <false>
            Smap2Hypix             = false
            Temporary              = false # <true>; <false>

			run = RUN(ChangeHydroModel, HydroLabθΨ⍰, Hypix, Infilt, IntergranularMixingPsd, RockCorection, Smap, Smap2Hypix, Temporary)
				
			# =============================================================
			#	   ROCK FRAGMENT OPTIONS
			# =============================================================
				# Rocks options
					CorectStoneRockWetability = false # <true> rocks are wettable; <false> rocks are mot wettable 
					RockInjectedIncluded⍰ = :InjectRock # <:InjectRock> rocks fragments are injected/forced into the fine soils; <Included> rocks are included in the bulk BulkDensity

				rockFragment = ROCKFRAGMENT(CorectStoneRockWetability, RockInjectedIncluded⍰)
				

			# =============================================================
			#		SMAP OPTIONS
			# =============================================================
				# Smap-Hydro options
					Plot_Kunsat = false  # <true> or <false>

				smap = SMAP(Plot_Kunsat)

			# =============================================================
			#		HYDRO OPTIONS
			# =============================================================
				# Hydraulic model
				HydroModel⍰      = :Kosugi # <:Kosugi>*; <:Vangenuchten>; <:BrooksCorey>; <:ClappHornberger>; <:VangenuchtenJules>
				# HydroModel_List = [:Kosugi; :Vangenuchten; :BrooksCorey; :ClappHornberger; :VangenuchtenJules]
				HydroModel_List = [:Kosugi]
					θsOpt⍰           = :Φ #  <:Opt> Optimize θs; <:Φ> derived from total porosity which requires some correction from param.hydro.Coeff_Φ_2_θs; <:FromData> θs is optimised with the feasible range derived directly from θobs(Ψ), assuming that we have θ(Ψ=0))
					θrOpt⍰           = :Opt # <:Opt> optimises; <:ParamPsd> derive from PSD and uses α1 and α1 from parameters in Param.jl; <:σ_2_θr>	
					σ_2_Ψm⍰          = :Constrained # <:Constrained> Ψm physical feasible range is computed from σ <:UniqueRelationship> Ψm is computed from σ; <:No> optimisation of σ & Ψm with no constraints

				# PLOTTING
					Plot_θΨ   = true
					Plot_σ_Ψm = false
					
				hydro = HYDRO(HydroModel⍰, HydroModel_List, θsOpt⍰, θrOpt⍰, σ_2_Ψm⍰, Plot_θΨ, Plot_σ_Ψm)


			# =============================================================
			#		PSD OPTIONS     
			# =============================================================
				# MODEL⍰
               Model⍰      = :IMP # <:IMP>* Intergranular Mixing Model⍰; <:Chang2019Model>
               OptimizePsd⍰ = :OptAllSoil # <:OptSingleSoil>; <:OptAllSoil>; or <:Run>
               Psd_2_θr⍰    = :Opt # <:Opt> optimises parameters α1 and α1; <:ParamPsd> uses α1 and α1 from parameters in Param.jl

				# For OptimizePsd⍰ = :single
					∑Psd_2_ξ1 = true  # <true> optimize ξ1; <false>
					
				# FITTING THE PSD FUNCTION TO A HYDRAULIC MODEL			
               HydroModel⍰ = :Kosugi # <:Kosugi>*; <:Vangenuchten>; <:BrooksCorey>; <:ClappHornberger>; <:VangenuchtenJules>
              
               θsOpt⍰      = :Φ #  <:Opt> Optimize θs; <:Φ> derived from total porosity which requires some correction from param.hydro.Coeff_Φ_2_θs; <:FromData> θs is optimised with the feasible range derived directly from θobs(Ψ), assuming that we have θ(Ψ=0))
               θrOpt⍰      = :Opt # <:Opt> optimises; <:ParamPsd> derive from PSD and uses α1 and α1 from parameters in Param.jl; <:σ_2_θr>
               σ_2_Ψm⍰     = :Constrained # <:Constrained> Ψm physical feasible range is computed from σ <:UniqueRelationship> Ψm is computed from σ; <:No> optimisation of σ & Ψm with no constraints

				# PLOTTING
					Plot_Psd_θΨ    = true # <true>  plot θΨ of PSD; <false>
					Plot_θr        = true #  <true>  plot θr data and model from Psd ; <false>	
					Plot_IMP_Model = true # <true> ; plot IMP model results for publication; <false>
					
				# TABLE
					Table_Psd_θΨ_θ = true # <true> derive θ values at prescribed Ψ

				psd = PSD(Model⍰, OptimizePsd⍰, Psd_2_θr⍰, ∑Psd_2_ξ1, HydroModel⍰, θsOpt⍰, θrOpt⍰, σ_2_Ψm⍰, Plot_Psd_θΨ, Plot_θr, Plot_IMP_Model, Table_Psd_θΨ_θ)


			# =============================================================
			#		INFILTRATION OPTIONS
			# =============================================================
				# Model⍰
               DataSingleDoubleRing⍰ = :Double	# <:Double> infiltration measured by double ring infiltrometer; <:Single> infiltration measured by single ring infiltrometer
               OptimizeRun⍰          = :Opt # <:Opt>* Optimise hydraulic parameters from infiltration data; <:Run> run inftration curves from known hydraulic parameters; <:RunOptKs>  run inftration curves from known hydraulic parameters but optimize Ks only
               Model⍰                = :QuasiExact 	# <:QuasiExact> physical approach; <:Best_Univ> statistical improved approach
               BestUniv_Continous    = true # <true> if the function is continous as recommended; <false> traditional
               SorptivitySplitModel⍰ = :Split # <:Split>; <:Split_η>
               SorptivityModel⍰      = :Parlange # <:Parlange> strong non-linear diffusivity;  <:Crank> constant diffusivity; <:Philip_Knight> dirac delta-function diffusivity; <:Brutsaert> moderate non-linear diffusivity,
				
				# Deriving hydraulic parameters from infiltration tests
					HydroModel⍰      = :Kosugi # <:Kosugi>*; <:Vangenuchten>; <:BrooksCorey>; <:ClappHornberger>; <:VangenuchtenJules>
				
					θsOpt⍰           = :Φ #  <:Opt> Optimize θs; <:Φ> derived from total porosity which requires some correction from param.hydro.Coeff_Φ_2_θs; <:FromData> θs is optimised with the feasible range derived directly from θobs(Ψ), assuming that we have θ(Ψ=0))
					θrOpt⍰           = :Opt # <:Opt> optimises; <:ParamPsd> derive from PSD and uses α1 and α1 from parameters in Param.jl; <:σ_2_θr>	
					σ_2_Ψm⍰          = :Constrained # <:Constrained> Ψm physical feasible range is computed from σ <:UniqueRelationship> Ψm is computed from σ; <:No> optimisation of σ & Ψm with no constraints
			
				# Plotting
					Plot_Sorptivity       = true # <true> or <false>	
					Plot_SeIni_Range      = true # <true> computes infiltration curves for different SeIn set in param.infilt.SeIni_Output <false> no outputs
					Plot_∑Infilt          = true # <true> plots cumulative infiltration curves for experimental and derived data <false> no plots
					Plot_θΨ               = true # <true>; <false>
					Plot_Sorptivity_SeIni = true # <true> computes sorptivity curves as a function of Se <false> no outputs

				infilt = INFILT(DataSingleDoubleRing⍰, OptimizeRun⍰, Model⍰, BestUniv_Continous, SorptivitySplitModel⍰, SorptivityModel⍰ , HydroModel⍰, θsOpt⍰, θrOpt⍰, σ_2_Ψm⍰, Plot_Sorptivity, Plot_SeIni_Range, Plot_∑Infilt, Plot_θΨ, Plot_Sorptivity_SeIni)


			# =============================================================
			#		HYPIX OPTIONS
			# =============================================================
				# Time step
					ClimateDataTimestep⍰ = "Daily" # <Hourly>; <Daily>

				# Modules used
					RainfallInterception = true
					Evaporation          = true
					RootWaterUptake      = true
					RootWaterUptakeComp  = true
					
				#S sink term 
					LookupTable_Lai            = true # <false> Lai=constant; <true> Lai varies per month
					LookUpTable_CropCoeficient = true # <false> CropCoeficient=constant; <true> CropCoeficient varies per month

				# Hydraulic model 
					HydroModel⍰ = :Kosugi # <:vanGenuchten>; <:Kosugi>

				# Richards equation
					BottomBoundary = :Free # not working <:Free>; <:Pressure>
					∂R∂Ψ_Numerical = false # perform the derivatives numerically <true>; <false>

				# Adaptive time step
					AdaptiveTimeStep   = :ΔΨ # <:ΔΨ>; <:Δθ>
					NormMin            = :Norm		#<:Norm>; <:Min>
					Flag_ReRun         = true # <true>; <false> Rerun after updating the ΔT
					Qbottom_Correction = true # <true> correction for the mass balance of the last cell
					# const NoConverge_Ψbest   = false # not working <true>; <false>* compute Q(Ψbest) when no convergence
				
				# Rainfall interception model
					Lai_2_SintMax = false # <true> derive Sint_Sat from LAI_2_SINTMAX; <false> derive from input file

				# Step wise optimization
					σ_2_Ψm⍰ = :No  # <:Constrained> Ψm physical feasible range is computed from σ <:UniqueRelationship> Ψm is computed from σ; <:No> optimisation of σ & Ψm with no constraints
					σ_2_θr = true # <true> derive θr from σ <false>
					θs_Opt = :No #  <:θs_Opt> θs is derived by multiplying a parameter to Max(θobs) for all profiles; <No>

				# Calibration data available
					Optimisation = true # <true>; <false>
					θobs         = true # <true>; <false>
					θobs_Average = true #<true> ; <false>determine if the observed θ is an average of different layers

					θobs_Hourly = true # θ data can be very large so we reduce the data to hourly
					Signature_Run = false

				# Table true
					Table = true 
					Table_Discretization  = true
					Table_Q               = true
					Table_RootWaterUptake = true
					Table_TimeSeries      = true
					Table_Ψ               = true
					Table_θ               = true
					Table_TimeSeriesDaily = true
					Tabule_θΨ             = true
					Table_Climate         = true
				
				# plot outputs
					Plot_Vegetation   = false
					Plot_θΨK          = false
					Plot_Interception = false
					Plot_Other        = false
					Plot_Sorptivity   = false
					Plot_Hypix        = true
						Plot_Climate      = true
						Plot_θ            = true
						Plot_Ψ            = true
						Plot_Flux         = true
						Plot_WaterBalance = true
						Plot_ΔT           = true

			hyPix = HYPIX(ClimateDataTimestep⍰, RainfallInterception, Evaporation, RootWaterUptake, RootWaterUptakeComp, LookupTable_Lai, LookUpTable_CropCoeficient, HydroModel⍰, BottomBoundary, ∂R∂Ψ_Numerical, AdaptiveTimeStep, NormMin, Flag_ReRun, Qbottom_Correction, Lai_2_SintMax, σ_2_Ψm⍰, σ_2_θr, θs_Opt, Optimisation, θobs,θobs_Average, θobs_Hourly, Signature_Run, Table, Table_Discretization, Table_Q, Table_RootWaterUptake, Table_TimeSeries, Table_Ψ, Table_θ, Table_TimeSeriesDaily, Tabule_θΨ, Table_Climate, Plot_Vegetation, Plot_θΨK, Plot_Interception, Plot_Other, Plot_Sorptivity, Plot_Hypix, Plot_Climate, Plot_θ, Plot_Ψ, Plot_Flux, Plot_WaterBalance, Plot_ΔT)


			# =============================================================
			#		GLOBAL OPTION
			# ===========================================================
				option = OPTION(data, dataFrom, hydro, hyPix, infilt, other, psd, rockFragment, run, smap)

		return option
		end  # function: OPTION

end # module option 
# end OPTION