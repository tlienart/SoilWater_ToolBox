# =============================================================
#		module option
# ===========================================================
module options

	using Configurations, TOML

	# What available data we have?
	@option struct DATA
		HydroParamPrecomputed::Bool
		Infilt::Bool
		Jules::Bool
		Kθ::Bool
		Nsdr::Bool
		Pedological⍰::String
		Psd::Bool
		RockWetability::Bool
		SimulationKosugiθΨK::Bool
		θΨ::Bool
		Φ⍰::String
	end # struct DATA

	# What model wanting to run
	@option mutable struct RUN
      ChangeHydroModel::Bool
      HydroLabθΨ⍰::String
      Hypix::Bool
      Infilt::Bool
      IntergranularMixingPsd ::Bool
		Jules::Bool
		KsModel::Bool
      RockCorection::Bool
		Smap::Bool
		Smap2Hypix::Bool
      Temporary::Bool
	end

	@option mutable struct OTHER
		DownloadPackage::Bool
		Ploting::Bool
		PlotVscode::Bool
	end
	@option mutable struct SMAP
		Nothings
	end
	@option mutable struct HYDRO
		HydroModel⍰::String
		HydroModel_List
		θsOpt⍰::String
		θrOpt⍰::String
		σ_2_Ψm⍰::String
		Plot_θΨ::Bool
	end

	@option mutable struct KSMODEL
		KsModel⍰::String
	end

	@option mutable struct PSD
		Model⍰::String
		OptimizePsd⍰::String
		Psd_2_θr⍰::String
		∑Psd_2_ξ1::Bool
		HydroModel⍰::String
		θsOpt⍰::String
		θrOpt⍰::String
		σ_2_Ψm⍰::String
		Plot_Psd_θΨ::Bool
		Plot_θr::Bool
		Plot_IMP_Model::Bool
		Table_Psd_θΨ_θ::Bool
	end
	@option mutable struct ROCKFRAGMENT
		CorectStoneRockWetability::Bool
		RockInjectedIncluded⍰::String
	end
	@option mutable struct INFILT
		DataSingleDoubleRing⍰::String  
		OptimizeRun⍰::String  
		Model⍰::String
		BestUniv_Continous::Bool                	  
		SorptivitySplitModel⍰::String  
		SorptivityModel⍰::String     
		HydroModel⍰::String       
		θsOpt⍰::String            
		θrOpt⍰::String            
		σ_2_Ψm⍰::String               
		Plot_Sorptivity::Bool        	 
		Plot_∑Infilt::Bool           
		Plot_θΨ::Bool                
	end
	@option mutable struct HYPIX
		ClimateDataTimestep⍰::String
		RainfallInterception::Bool
		Evaporation::Bool
		RootWaterUptake::Bool
		RootWaterUptakeComp::Bool
		LookupTable_Lai::Bool
		LookUpTable_CropCoeficient::Bool
		HydroModel⍰::String
		BottomBoundary⍰
		∂R∂Ψ_Numerical::Bool
		AdaptiveTimeStep⍰::String
		NormMin⍰::String
		Flag_ReRun::Bool
		Qbottom_Correction::Bool
		Lai_2_SintMax::Bool
		σ_2_Ψm⍰::String
		σ_2_θr::Bool
		θs_Opt⍰::String
		Optimisation::Bool
		θobs::Bool
		θobs_Average::Bool
		θobs_Hourly::Bool
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
		@option mutable struct OPTION
			data::DATA
			hydro::HYDRO
			hyPix::HYPIX
			infilt::INFILT
			other::OTHER
			psd::PSD
			rockFragment::ROCKFRAGMENT
			run::RUN
			smap::SMAP
			ksModel::KSMODEL
		end # struct OPTION
	
	#__________________________________________________________________
	#..................................................................

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPTION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OPTIONS(Path_Data, SiteName)
			
			Path = Path_Data * "/ParamOptionPath/" * SiteName * "_Option.toml"

			return option = Configurations.from_toml(OPTION, Path)
		end  # function: OPTION

end # module option 
# end OPTION