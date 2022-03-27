# =============================================================
#		module optionHypix hypix
# =============================================================

"""
   OPTION_HYPIX 
automatically puts the values of options from toml file into mutuable structures optionHypix
"""
module optionsHypix

	using Configurations, TOML

	@option mutable struct OPT
		σ_2_Ψm⍰::String
		σ_2_θr::Bool
		θs_Opt⍰::String
		Optimisation::Bool
		HydroVegParamReadFromOutput::Bool
	end

	@option mutable struct OPTIONHYPIX
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
		∂R∂Ψ_NumericalAuto::Bool
		# IterReduceOverShoting::Bool
		# DynamicNewtonRaphsonStep::Bool
		AdaptiveTimeStep⍰::String
		# Flag_ReRun::Bool
		Lai_2_SintMax::Bool

		θobs::Bool
		θavr_RootZone::Bool
		θobs_Reduced::Bool
		
		Table::Bool
		Table_Discretization::Bool
		Table_Q::Bool
		Table_Ψ::Bool
		Table_θ::Bool
		Table_TimeSerie::Bool
		Tabule_θΨ::Bool
		Ploting::Bool
		Plot_θprofile::Bool
		Plot_Vegetation::Bool
		Plot_θΨK::Bool
		Plot_Interception::Bool
		Plot_Other::Bool
		Plot_Sorptivity::Bool
		Plot_Hypix::Bool
		Plot_Climate::Bool
		Plot_θ::Bool
		Plot_Flux::Bool
		opt::OPT
	end


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPTION_HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OPTION_HYPIX(PathOptionHypix::String)
			return Configurations.from_toml(OPTIONHYPIX, PathOptionHypix)
		end  # function: OPTION_HYPIX

end # module optionHypix
#.................................................................. 