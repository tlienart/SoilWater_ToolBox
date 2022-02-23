# =============================================================
#		MODULE: pathOutputHypix
# =============================================================
module pathsHypix

	using Configurations

	@option mutable struct PATHYPIX
		Input_OfStep::String

		Table_DailyClimate::String
		Table_Discretisation::String
		Table_Hydro::String
		Table_KΨ::String
		Table_Performance::String
		Table_Q::String
		Table_Signature::String
		Table_TimeSerie::String
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
		Plot_Se_Z::String
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
				pathOutputHypix = Configurations.from_toml(PATHYPIX, PathPathHypix)

			# Paths Output
				Path_OutputHypix = Path_Hypix * "\\data\\OUTPUT\\Hypix\\" * ProjectHypix * "\\" * SiteName₀ *  "\\" 
				
			# HYPIX OUTPUT TABLE
				Path_Hypix_Table = Path_OutputHypix *"/Table/" 

					mkpath(Path_Hypix_Table) #Make Folder if not exist

					Path_Hypix_Table = Path_Hypix_Table * SiteName₀ * "_"

					pathOutputHypix.Table_DailyClimate    = Path_Hypix_Table  *  SiteName₀ * "_"* pathOutputHypix.Table_DailyClimate
					pathOutputHypix.Table_Discretisation  = Path_Hypix_Table  *  SiteName₀ * "_" *pathOutputHypix.Table_Discretisation
					pathOutputHypix.Table_Hydro           = Path_Hypix_Table  *  SiteName₀ * "_" *pathOutputHypix.Table_Hydro
					pathOutputHypix.Table_KΨ              = Path_Hypix_Table  *  SiteName₀ * "_"* pathOutputHypix.Table_KΨ
					pathOutputHypix.Table_Performance     = Path_Hypix_Table  *  SiteName₀ * "_"* pathOutputHypix.Table_Performance
					pathOutputHypix.Table_Q               = Path_Hypix_Table  *  SiteName₀ * "_"* pathOutputHypix.Table_Q
					pathOutputHypix.Table_Signature       = Path_Hypix_Table  *  SiteName₀ * "_"* pathOutputHypix.Table_Signature
					pathOutputHypix.Table_TimeSerie       = Path_Hypix_Table  *  SiteName₀ * "_"* pathOutputHypix.Table_TimeSerie
					pathOutputHypix.Table_TimeSerie_Daily = Path_Hypix_Table  *  SiteName₀ * "_"* pathOutputHypix.Table_TimeSerie_Daily
					pathOutputHypix.Table_Veg             = Path_Hypix_Table  *  SiteName₀ * "_"* pathOutputHypix.Table_Veg
					pathOutputHypix.Table_θ               = Path_Hypix_Table  *  SiteName₀ * "_"* pathOutputHypix.Table_θ
					pathOutputHypix.Table_θΨ              = Path_Hypix_Table  *  SiteName₀ * "_"* pathOutputHypix.Table_θΨ
					pathOutputHypix.Table_Ψ               = Path_Hypix_Table  *  SiteName₀ * "_"* pathOutputHypix.Table_Ψ

				# HYPIX PLOT
					Path_Hypix_Plot = Path_OutputHypix * "/Plots/" 	
						mkpath(Path_Hypix_Plot)

						pathOutputHypix.Plot_HypixTime            = Path_Hypix_Plot * SiteName₀  * "_" * pathOutputHypix.Plot_HypixTime
						pathOutputHypix.Plot_θprofile             = Path_Hypix_Plot * SiteName₀  * "_" * pathOutputHypix.Plot_θprofile
						pathOutputHypix.Plot_Hypix_θΨK            = Path_Hypix_Plot * SiteName₀  * "_" * pathOutputHypix.Plot_Hypix_θΨK
						pathOutputHypix.Plot_RainfallInterception = Path_Hypix_Plot * SiteName₀  * "_" * pathOutputHypix.Plot_RainfallInterception
						pathOutputHypix.Plot_Se_Time              = Path_Hypix_Plot * SiteName₀  * "_" * pathOutputHypix.Plot_Se_Time
						pathOutputHypix.Plot_Se_Z                 = Path_Hypix_Plot * SiteName₀  * "_" * pathOutputHypix.Plot_Se_Z
						pathOutputHypix.Plot_Sorptivity           = Path_Hypix_Plot * SiteName₀  * "_" * pathOutputHypix.Plot_Sorptivity
						pathOutputHypix.Vegetation                = Path_Hypix_Plot * SiteName₀  * "_" * pathOutputHypix.Vegetation

				# HYPIX OTHERS
					Path_Hypix_Other =  Path_Hypix * "\\data\\OUTPUT\\Hypix\\" * ProjectHypix * "\\" * "OTHER" *  "\\"  
						mkpath(Path_Hypix_Other)

						pathOutputHypix.Plot_OfStep           = Path_Hypix_Other
						pathOutputHypix.Plot_θ∂θ∂Ψ            = Path_Hypix_Other* pathOutputHypix.Plot_θ∂θ∂Ψ
						pathOutputHypix.Plot_Ψmin_Ψmax        = Path_Hypix_Other* pathOutputHypix.Plot_Ψmin_Ψmax
						pathOutputHypix.Plot_σ2θr             = Path_Hypix_Other * pathOutputHypix.Plot_σ2θr
						# pathOutputHypix.Plot_θΨ_Δθ            = Path_Hypix_Other * pathOutputHypix.Plot_θΨ_Δθ
						pathOutputHypix.Plot_Se_Ψ_Constrained = Path_Hypix_Other * pathOutputHypix.Plot_Se_Ψ_Constrained

		return pathsHypix
		end # function PATH_HYPIX()	

end  # module pathsHypix
# ............................................................