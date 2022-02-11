# =============================================================
#		module option hypix
# =============================================================

"""
   OPTION_HYPIX 
automatically puts the values of options from toml file into mutuable structures optionHypix
"""
module optionHypix

	using Configurations, TOML

	@option mutable struct OTHER
      Ploting::Bool
      PlotVscode::Bool
	end

	@option mutable struct HYPIX
		RainfallInterception::Bool
		Evaporation::Bool
		RootWaterUptake::Bool
		RootWaterUptakeComp::Bool
		Ponding::Bool
		LookupTable_Lai::Bool
		LookUpTable_CropCoeficient::Bool
		Discretisation_File_Auto⍰::String
		HydrostaticEquilibrium::Bool
		HydroModel⍰::String
		TopBoundary⍰::String
		BottomBoundary⍰::String
		∂R∂Ψ_Numerical::Bool
		IterReduceOverShoting::Bool
		DynamicNewtonRaphsonStep::Bool
		ZhaWetingDrySoil::Bool
		AdaptiveTimeStep⍰::String
		NormMin⍰::String
		Flag_ReRun::Bool
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
		Plot_θprofile::Bool
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

   @option mutable struct OPTIONHYPIX
		other::OTHER
      hypix::HYPIX
   end # struct OPTION
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPTION_HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OPTION_HYPIX(iSim, PathData_Hypix, SiteName_Hypix, Soilname)

			# Specific path for which the options are only true for the Soilname in question
			Path = PathData_Hypix * "/" * Soilname[iSim] *  "/ParamOptionPath/" * SiteName_Hypix * "_OptionHypix.toml" 

			if !isfile(Path)
				# Global path for which the options are true for all Soilname
				Path = PathData_Hypix * "/" *  Soilname[iSim] * "/ParamOptionPath/" * Soilname[iSim] * "_OptionHypix.toml"

			elseif !isfile(Path)
				error("Cannot find $Path")
			end

		return Configurations.from_toml(OPTIONHYPIX, Path)
		end  # function: OPTION_HYPIX

end # module optionHypix
#.................................................................. 
