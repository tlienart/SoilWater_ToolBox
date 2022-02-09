# =============================================================
#		MODULE: pathHypix
# =============================================================
module pathHypix

	using Configurations

	@option struct OPTIONS 
		ModelName::String
		Select::String
	end # struct OPTION

	@option mutable struct PATHYPIX
		Climate::String
		Dates::String
		Discretization::String
		DiscretizationAuto::String
		HyPix_HydroParam::String
		HyPixParamOpt::String

		HyPix_VegParam::String
		IdSelect::String 
		Input_OfStep::String
		JulesMetadata::String
		IdName_Hypix::String
		obsTheta::String 

		LookUpTable_CropCoeficient::String
		LookUpTable_Lai::String

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
		Table_θaverage::String
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

	@option mutable struct PATHHYPIX
		hyPix::PATHYPIX
		option::OPTIONS
	end

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PATH_HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PATH_HYPIX(iSim, PathData_Hypix, SiteName_Hypix; Soilname=["Temporary"])

			# Specific path for which the options are only true for the Soilname in question
			Path = PathData_Hypix * "/" * Soilname[iSim] *  "/ParamOptionPath/" * SiteName_Hypix * "_PathHypix.toml" 

			if !isfile(Path)
				# Global path for which the options are true for all Soilname
				Path = PathData_Hypix * "/" *  Soilname[iSim] * "/ParamOptionPath/" * Soilname[iSim] * "_PathHypix.toml"

			elseif !isfile(Path)
				error("Cannot find $Path")
			end

		# Reading toml -> path structure
			path = Configurations.from_toml(PATHHYPIX, Path)

		# Usefull path
			Path_Home₀ = @__DIR__

		# perform cs..
			Path_Home = dirname(Path_Home₀)
			Path_Home = Path_Home * "/data/"
			
		
		# =============================================================
		#		HYPIX MODEL
		# =============================================================
			IdName_Hypix = Soilname[iSim]
		
		# HYPIX INPUT LEVEL 1 ===
         Path_Hypix₁         = PathData_Hypix * "/"
         path.hyPix.IdSelect      = Path_Hypix₁ * SiteName_Hypix * "_" * path.hyPix.IdSelect
         path.hyPix.JulesMetadata = Path_Hypix₁ * path.hyPix.JulesMetadata

		# HYPIX INPUT LEVEL 2 ===
		Path_Hypix₂            = PathData_Hypix * "/" * IdName_Hypix * "/" * IdName_Hypix * "_"

         path.hyPix.Climate            = Path_Hypix₂ * path.hyPix.Climate
         path.hyPix.Dates              = Path_Hypix₂ * path.hyPix.Dates
         path.hyPix.Discretization     = Path_Hypix₂ * path.hyPix.Discretization
         path.hyPix.DiscretizationAuto = Path_Hypix₂ * path.hyPix.DiscretizationAuto
         path.hyPix.HyPix_HydroParam   = Path_Hypix₂ * path.hyPix.HyPix_HydroParam
         path.hyPix.HyPix_VegParam     = Path_Hypix₂ * path.hyPix.HyPix_VegParam
         path.hyPix.HyPixParamOpt      = Path_Hypix₂ * path.hyPix.HyPixParamOpt
         path.hyPix.obsTheta           = Path_Hypix₂ * path.hyPix.obsTheta

			path.hyPix.Input_OfStep     = PathData_Hypix * "/RESULTS/"

		# HYPIX LOOKUPTABLE ===
      Path_Hypix_LookUpTable                 = Path_Home * "/INPUT/Data_Hypix/LookUpTable/"
					
         path.hyPix.LookUpTable_CropCoeficient = Path_Hypix_LookUpTable * path.hyPix.LookUpTable_CropCoeficient
         path.hyPix.LookUpTable_Lai                       = Path_Hypix_LookUpTable * path.hyPix.LookUpTable_Lai

			# HYPIX OUTPUT TABLE
			Path_Hypix_Table = Path_Home * "/OUTPUT/Hypix/" * SiteName_Hypix * "/" * IdName_Hypix *"/Table/" 				
				mkpath(Path_Hypix_Table) #Make Folder if not exist

				Path_Hypix_Table = Path_Hypix_Table * SiteName_Hypix * "_"

            path.hyPix.Table_DailyClimate    = Path_Hypix_Table  *  IdName_Hypix * "_"* path.hyPix.Table_DailyClimate
            path.hyPix.Table_Discretisation  = Path_Hypix_Table  *  IdName_Hypix * "_" *path.hyPix.Table_Discretisation
            path.hyPix.Table_Hydro           = Path_Hypix_Table  *  IdName_Hypix * "_" *path.hyPix.Table_Hydro
            path.hyPix.Table_KΨ              = Path_Hypix_Table  *  IdName_Hypix * "_"* path.hyPix.Table_KΨ
            path.hyPix.Table_Performance     = Path_Hypix_Table  *  IdName_Hypix * "_"* path.hyPix.Table_Performance
            path.hyPix.Table_Q               = Path_Hypix_Table  *  IdName_Hypix * "_"* path.hyPix.Table_Q
            path.hyPix.Table_Signature       = Path_Hypix_Table  *  IdName_Hypix * "_"* path.hyPix.Table_Signature
            path.hyPix.Table_TimeSerie       = Path_Hypix_Table  *  IdName_Hypix * "_"* path.hyPix.Table_TimeSerie
            path.hyPix.Table_TimeSerie_Daily = Path_Hypix_Table  *  IdName_Hypix * "_"* path.hyPix.Table_TimeSerie_Daily
            path.hyPix.Table_Veg             = Path_Hypix_Table  *  IdName_Hypix * "_"* path.hyPix.Table_Veg
            path.hyPix.Table_θ               = Path_Hypix_Table  *  IdName_Hypix * "_"* path.hyPix.Table_θ
            path.hyPix.Table_θΨ              = Path_Hypix_Table  *  IdName_Hypix * "_"* path.hyPix.Table_θΨ
            path.hyPix.Table_Ψ               = Path_Hypix_Table  *  IdName_Hypix * "_"* path.hyPix.Table_Ψ
				
			Path_Hypix_Table_θaverage = Path_Home * "/OUTPUT/Hypix/" * SiteName_Hypix * "/SoilMoistureSim/" 				
				mkpath(Path_Hypix_Table_θaverage) #Make Folder if not exist
				path.hyPix.Table_θaverage        = Path_Hypix_Table_θaverage *  IdName_Hypix * "_"* path.hyPix.Table_θaverage

			# HYPIX PLOT CORE
			Path_Hypix_Plot = Path_Home * "/OUTPUT/Hypix/" * SiteName_Hypix * "/" * IdName_Hypix *"/Plots/" 	
				mkpath(Path_Hypix_Plot)

            path.hyPix.Plot_HypixTime            = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_HypixTime
            path.hyPix.Plot_θprofile             = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_θprofile
            path.hyPix.Plot_Hypix_θΨK            = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Hypix_θΨK
            path.hyPix.Plot_RainfallInterception = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_RainfallInterception
            path.hyPix.Plot_Se_Time              = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Se_Time
            path.hyPix.Plot_Se_Z                 = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Se_Z
            path.hyPix.Plot_Sorptivity           = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Sorptivity
            path.hyPix.Vegetation                = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Vegetation

			# HYPIX PLOT OTHERS: RESULTS
			Path_Hypix_Result= Path_Home * "/OUTPUT/Hypix/RESULTS/"
				mkpath(Path_Hypix_Result)

				path.hyPix.Plot_OfStep           = Path_Hypix_Result
				path.hyPix.Plot_θ∂θ∂Ψ            = Path_Hypix_Result* path.hyPix.Plot_θ∂θ∂Ψ
				path.hyPix.Plot_Ψmin_Ψmax        = Path_Hypix_Result* path.hyPix.Plot_Ψmin_Ψmax
				path.hyPix.Plot_σ2θr             = Path_Hypix_Result * path.hyPix.Plot_σ2θr
				path.hyPix.Plot_θΨ_Δθ            = Path_Hypix_Result * path.hyPix.Plot_θΨ_Δθ
				path.hyPix.Plot_Se_Ψ_Constrained = Path_Hypix_Result * path.hyPix.Plot_Se_Ψ_Constrained

	return path
	end # function PATHS			
end  # module path
# ............................................................