# =============================================================
#		module option
# ===========================================================
module options

	using Configurations, TOML

	# What available data we have?
	@option struct DATA
		HydroParam::Bool
		Infilt::Bool
		Jules::Bool
		Kθ::Bool
		Pedological⍰
		Psd::Bool
		RockWetability::Bool
		SimulationKosugiθΨK::Bool
		θΨ::Bool
		Φ⍰
	end # struct DATA

	# What model wanting to run
	@option mutable struct RUN
      ChangeHydroModel::Bool
      HydroLabθΨ⍰
      Hypix::Bool
      Infilt::Bool
      IntergranularMixingPsd ::Bool
		Jules::Bool
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
		HydroModel⍰
		HydroModel_List
		θsOpt⍰
		θrOpt⍰
		σ_2_Ψm⍰
		Plot_θΨ::Bool
	end
	@option mutable struct PSD
		Model⍰
		OptimizePsd⍰
		Psd_2_θr⍰
		∑Psd_2_ξ1::Bool
		HydroModel⍰
		θsOpt⍰
		θrOpt⍰
		σ_2_Ψm⍰
		Plot_Psd_θΨ::Bool
		Plot_θr::Bool
		Plot_IMP_Model::Bool
		Table_Psd_θΨ_θ::Bool
	end
	@option mutable struct ROCKFRAGMENT
		CorectStoneRockWetability::Bool
		RockInjectedIncluded⍰
	end
	@option mutable struct INFILT
		DataSingleDoubleRing⍰  
		OptimizeRun⍰  
		Model⍰
		BestUniv_Continous::Bool                	  
		SorptivitySplitModel⍰  
		SorptivityModel⍰     
		HydroModel⍰       
		θsOpt⍰            
		θrOpt⍰            
		σ_2_Ψm⍰               
		Plot_Sorptivity::Bool        	 
		Plot_∑Infilt::Bool           
		Plot_θΨ::Bool                
	end
	@option mutable struct HYPIX
		ClimateDataTimestep⍰
		RainfallInterception::Bool
		Evaporation::Bool
		RootWaterUptake::Bool
		RootWaterUptakeComp::Bool
		LookupTable_Lai::Bool
		LookUpTable_CropCoeficient::Bool
		HydroModel⍰
		BottomBoundary⍰
		∂R∂Ψ_Numerical::Bool
		AdaptiveTimeStep⍰
		NormMin⍰
		Flag_ReRun::Bool
		Qbottom_Correction::Bool
		Lai_2_SintMax::Bool
		σ_2_Ψm⍰
		σ_2_θr::Bool
		θs_Opt⍰
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
		end # struct OPTION
	
	#__________________________________________________________________
	#..................................................................

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPTION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OPTIONS()

			Path= "D:\\Main\\MODELS\\SoilWater-ToolBox2\\src\\Option.toml"

			TomlOption = TOML.parsefile(Path)



			option = Configurations.from_dict(OPTION, TomlOption)

			FieldsInStruct=fieldnames(typeof(option));
			for i=1:length(FieldsInStruct)
				#Check field i
				Value=getfield(option, FieldsInStruct[i])
				# Value=getfield(Value, FieldsInStruct[i])
				println(Value)
				# Value=Value.+1;
				# setfield!(option,FieldsInStruct[i],Value)
			end

		return option
		end  # function: OPTION

end # module option 
# end OPTION

options.OPTIONS()