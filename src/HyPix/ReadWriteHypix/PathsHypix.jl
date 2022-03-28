# =============================================================
#		MODULE: pathOutputHypix
# =============================================================
module pathsHypix

	using Configurations

	@option mutable struct PATHYPIX
		ScenarioName::String

		Input_OfStep::String

		Table_Discretisation::String
		Table_Hydro::String
		Table_KΨ::String
		Table_Performance::String
		Table_Q::String
		Table_TimeSerie_Daily::String
		Table_Veg::String
		Table_Ψ::String
		Table_θ::String
		Table_θΨ::String

		Plot_Hypix_θΨK::String
		Plot_HypixTime::String
		Plot_θprofile::String
		Plot_RainfallInterception::String
		Plot_Se_Time::String
		Plot_Sorptivity::String
		Vegetation::String

		Plot_OfStep::String
		Plot_Se_Ψ_Constrained::String
		Plot_θΨ_Δθ::String
		Plot_σ2θr::String
		Plot_Ψmin_Ψmax::String
		Plot_θ∂θ∂Ψ::String
	end

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PATH_HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PATH_HYPIX(Path_Hypix::String, PathPathHypix::String, ProjectHypix::String, SiteName₀::SubString{String})
				
			# READING toml -> path structure
				pathHypix = Configurations.from_toml(PATHYPIX, PathPathHypix)

			# Paths Output
				Path_OutputHypix = Path_Hypix * "\\data\\OUTPUT\\Hypix\\" * ProjectHypix * "\\" * SiteName₀ *  "\\" 
				
			# HYPIX OUTPUT TABLE
				Path_Hypix_Table = Path_OutputHypix *"/Table/" 

					mkpath(Path_Hypix_Table) #Make Folder if not exist

					Path_Hypix_Table = Path_Hypix_Table * SiteName₀ * "_"

					pathHypix.Table_Discretisation  = Path_Hypix_Table  *  pathHypix.ScenarioName * "_" *pathHypix.Table_Discretisation
					pathHypix.Table_Hydro           = Path_Hypix_Table  *  pathHypix.ScenarioName * "_" *pathHypix.Table_Hydro
					pathHypix.Table_KΨ              = Path_Hypix_Table  *  pathHypix.ScenarioName * "_"* pathHypix.Table_KΨ
					pathHypix.Table_Performance     = Path_Hypix_Table  *  pathHypix.ScenarioName * "_"* pathHypix.Table_Performance
					pathHypix.Table_Q               = Path_Hypix_Table  *  pathHypix.ScenarioName * "_"* pathHypix.Table_Q
					pathHypix.Table_TimeSerie_Daily = Path_Hypix_Table  *  pathHypix.ScenarioName * "_"* pathHypix.Table_TimeSerie_Daily
					pathHypix.Table_Veg             = Path_Hypix_Table  *  pathHypix.ScenarioName * "_"* pathHypix.Table_Veg
					pathHypix.Table_θ               = Path_Hypix_Table  *  pathHypix.ScenarioName * "_"* pathHypix.Table_θ
					pathHypix.Table_θΨ              = Path_Hypix_Table  *  pathHypix.ScenarioName * "_"* pathHypix.Table_θΨ
					pathHypix.Table_Ψ               = Path_Hypix_Table  *  pathHypix.ScenarioName * "_"* pathHypix.Table_Ψ

				# HYPIX PLOT
					Path_Hypix_Plot = Path_OutputHypix * "/Plots/" 	
						mkpath(Path_Hypix_Plot)

						Path_Hypix_Plot = Path_Hypix_Plot * SiteName₀ * "_"

						pathHypix.Plot_HypixTime            = Path_Hypix_Plot * pathHypix.ScenarioName  * "_" * pathHypix.Plot_HypixTime
						pathHypix.Plot_θprofile             = Path_Hypix_Plot * pathHypix.ScenarioName  * "_" * pathHypix.Plot_θprofile
						pathHypix.Plot_Hypix_θΨK            = Path_Hypix_Plot * pathHypix.ScenarioName  * "_" * pathHypix.Plot_Hypix_θΨK
						pathHypix.Plot_RainfallInterception = Path_Hypix_Plot * pathHypix.ScenarioName  * "_" * pathHypix.Plot_RainfallInterception
						pathHypix.Plot_Se_Time              = Path_Hypix_Plot * pathHypix.ScenarioName  * "_" * pathHypix.Plot_Se_Time
						pathHypix.Plot_Sorptivity           = Path_Hypix_Plot * pathHypix.ScenarioName  * "_" * pathHypix.Plot_Sorptivity
						pathHypix.Vegetation                = Path_Hypix_Plot * pathHypix.ScenarioName  * "_" * pathHypix.Vegetation

				# HYPIX OTHERS
					Path_Hypix_Other =  Path_Hypix * "\\data\\OUTPUT\\Hypix\\" * ProjectHypix * "\\" * "OTHER" *  "\\"  
						mkpath(Path_Hypix_Other)

						pathHypix.Plot_OfStep           = Path_Hypix_Other
						pathHypix.Plot_θ∂θ∂Ψ            = Path_Hypix_Other* pathHypix.Plot_θ∂θ∂Ψ
						pathHypix.Plot_Ψmin_Ψmax        = Path_Hypix_Other* pathHypix.Plot_Ψmin_Ψmax
						pathHypix.Plot_σ2θr             = Path_Hypix_Other * pathHypix.Plot_σ2θr
						# pathHypix.Plot_θΨ_Δθ            = Path_Hypix_Other * pathHypix.Plot_θΨ_Δθ
						pathHypix.Plot_Se_Ψ_Constrained = Path_Hypix_Other * pathHypix.Plot_Se_Ψ_Constrained

		return pathHypix
		end # function PATH_HYPIX()	

end  # module pathsHypix
# ............................................................