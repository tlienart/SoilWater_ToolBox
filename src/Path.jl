# =============================================================
#		MODULE: path
# =============================================================
module paths

	using Configurations

	@option struct OPTIONS 
		ModelName::String
		Select::String
	end # struct OPTION

	@option mutable struct INPUT_SMAP
		LookupTable_RockWetability::String
		Smap::String
	end

	@option mutable struct INPUT_SOILWATER
		BulkDensity::String
		ConvertModel::String
		HydroParam_Infilt::String
		HydroParam_ThetaH::String 
		IdSelect::String
		Infiltration::String
		Infiltration_Param::String
		Kunsat::String
		Psd::String
		Pedological⍰::String
		Φ::String
		Ψθ::String
	end # struct INPUT_SOILWATER

	@option mutable struct INPUT_TEMPORARY
		σ_ψM_Scenario::String
	end

	@option mutable struct TABLE_SMAP
		Table_Smap::String
		Table_θΨK::String
	end # struct INPUT_TABLE_SMAP	

	@option mutable struct TABLE_SOILWATER
		FileSoilHydro_Table₁::String
		Table_HydroInfilt::String
		Table_Infilt::String
		Table_KΨ::String
		Table_Psd_θΨ_θ::String
		Table_Psd::String
		Table_θΨ_Psd::String
		Table_θΨK::String
		TableComplete_KΨ::String
		TableComplete_θΨ::String
	end # struct INPUT_SOILWATER

	@option mutable struct PLOT_SOILWATER	
		Plot_∑infilt_Opt::String
		Plot_∑infilt_θΨ::String
		Plot_IMP_model::String
		Plot_Psd_θr::String
		Plot_Psd_θΨ::String				
		Plot_θΨK::String
		Plot_σΨm::String
	end # PLOT_SOILWATER

	@option mutable struct SMAP_2_HYPIX
		Path_Smap2Hypix::String
	end

	# _______________________ START: hypix _______________________ 
	@option mutable struct PATHYPIXS
		Climate::String
		Dates::String
		Discretization::String
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
# ------------------------END: hypix---------------------------  
@option mutable struct PATHS
		hyPix::PATHYPIXS
		inputSmap::INPUT_SMAP
		inputSoilwater::INPUT_SOILWATER
		inputTemporary::INPUT_TEMPORARY
		option::OPTIONS
		plotSoilwater::PLOT_SOILWATER
		smap2Hypix::SMAP_2_HYPIX
		tableSmap::TABLE_SMAP
		tableSoilwater::TABLE_SOILWATER
	end
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PATHS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PATH(iSim, opt, Path_Data, SiteName; Soilname=["Temporary"])

		PathToml = Path_Data * "/ParamOptionPath/" * SiteName * "_Path.toml"

		path = Configurations.from_toml(PATHS, PathToml)

		Path_Home₀ = @__DIR__

		# perform cs..
			Path_Home = dirname(Path_Home₀)

		# Change path name to /data/Private/
			Path_Home = Path_Home * "/data/"
			
		# =============================================================
		#		OPTIONS
		# =============================================================
			FileDataSoilhydro_Input = Path_Home * "INPUT/Data_SoilWater/" * SiteName * "/" * SiteName * "_"

			FileSoilHydro_Table₁ = Path_Home * "/OUTPUT/SoilWater/" * SiteName * "/Table/" 
				mkpath(FileSoilHydro_Table₁) 

		# =============================================================
		#		INPUT_SOILWATER
		# =============================================================
         path.inputSoilwater.BulkDensity        = FileDataSoilhydro_Input * path.inputSoilwater.BulkDensity
         path.inputSoilwater.ConvertModel       = FileDataSoilhydro_Input * path.inputSoilwater.ConvertModel
         path.inputSoilwater.HydroParam_Infilt  = FileDataSoilhydro_Input * path.inputSoilwater.HydroParam_Infilt
         path.inputSoilwater.HydroParam_ThetaH  = FileDataSoilhydro_Input * path.inputSoilwater.HydroParam_ThetaH
         path.inputSoilwater.IdSelect           = FileDataSoilhydro_Input * path.inputSoilwater.IdSelect
         path.inputSoilwater.Infiltration       = FileDataSoilhydro_Input * path.inputSoilwater.Infiltration
         path.inputSoilwater.Infiltration_Param = FileDataSoilhydro_Input * path.inputSoilwater.Infiltration_Param
         path.inputSoilwater.Kunsat             = FileDataSoilhydro_Input * path.inputSoilwater.Kunsat
         path.inputSoilwater.Psd                = FileDataSoilhydro_Input * path.inputSoilwater.Psd
         path.inputSoilwater.Pedological⍰      = FileDataSoilhydro_Input * path.inputSoilwater.Pedological⍰
         path.inputSoilwater.Φ                  = FileDataSoilhydro_Input * path.inputSoilwater.Φ
         path.inputSoilwater.Ψθ                 = FileDataSoilhydro_Input * path.inputSoilwater.Ψθ
		
		# =============================================================
		#		INPUT_SMAP
		#		None core
		# =============================================================
         path.inputSmap.LookupTable_RockWetability = FileDataSoilhydro_Input * path.inputSmap.LookupTable_RockWetability
         path.inputSmap.Smap                       = FileDataSoilhydro_Input * path.inputSmap.Smap

		# =============================================================
		#		INPUT_TEMPORARY
		#		None core
		# =============================================================
         path.inputTemporary.σ_ψM_Scenario = FileDataSoilhydro_Input * path.inputTemporary.σ_ψM_Scenario

		# =============================================================
		#		TABLE_SOILWATER
		# =============================================================
 
      path.tableSoilwater.FileSoilHydro_Table₁ = FileSoilHydro_Table₁ * SiteName

      path.tableSoilwater.Table_HydroInfilt    = FileSoilHydro_Table₁ * string(opt.infilt.Model⍰) * "_" *  path.option.ModelName  *  "_" * path.tableSoilwater.Table_HydroInfilt
      path.tableSoilwater.Table_Infilt         = FileSoilHydro_Table₁ * string(opt.infilt.Model⍰) *  "_" *  path.option.ModelName  *  "_" *  path.tableSoilwater.Table_Infilt
      path.tableSoilwater.Table_KΨ             = FileSoilHydro_Table₁ * "_"  *  path.option.ModelName * "_" * path.tableSoilwater.Table_KΨ
      path.tableSoilwater.Table_Psd            = FileSoilHydro_Table₁ * string(opt.psd.Model⍰) *  "_" * path.option.ModelName * "_" * path.tableSoilwater.Table_Psd
      path.tableSoilwater.Table_Psd_θΨ_θ       = FileSoilHydro_Table₁ * string(opt.psd.HydroModel⍰) *  "_" * path.option.ModelName * "_" *  path.tableSoilwater.Table_Psd_θΨ_θ
      path.tableSoilwater.Table_θΨ_Psd         = FileSoilHydro_Table₁ * string(opt.psd.HydroModel⍰) *  "_" * string(opt.hydro.σ_2_Ψm⍰) *  "_" * path.option.ModelName * "_" * path.tableSoilwater.Table_θΨ_Psd
      path.tableSoilwater.Table_θΨK            = FileSoilHydro_Table₁ *  "_" * string(opt.hydro.HydroModel⍰) * "_" * path.tableSoilwater.Table_θΨK
      path.tableSoilwater.TableComplete_θΨ     = FileSoilHydro_Table₁ *   "_" *  path.option.ModelName * "_" * path.tableSoilwater.TableComplete_θΨ
      path.tableSoilwater.TableComplete_KΨ     = FileSoilHydro_Table₁ *   "_" *  path.option.ModelName * "_" * path.tableSoilwater.TableComplete_KΨ

		# =============================================================
		#		TABLE_SMAP
		# =============================================================
         path.tableSmap.Table_θΨK  = FileSoilHydro_Table₁ * "_" * string(opt.hydro.HydroModel⍰) * "_" *  path.tableSmap.Table_θΨK
         path.tableSmap.Table_Smap = FileSoilHydro_Table₁ * path.tableSmap.Table_Smap

		# =============================================================
		#		PATH SMAP_2_HYPIX
		# =============================================================
			path.smap2Hypix.Path_Smap2Hypix = Path_Home *"OUTPUT/Smap2Hypix"
				mkpath(path.smap2Hypix.Path_Smap2Hypix)

	
		# =============================================================
		#		PLOT SOILWATER
		# =============================================================
			FileSoilHydro_Plot = Path_Home * "/OUTPUT/SoilWater/" * SiteName * "/Plots/"

			Plot_θΨK = FileSoilHydro_Plot * "/Lab/" 
				mkpath(Plot_θΨK)
				path.plotSoilwater.Plot_θΨK = Plot_θΨK * SiteName * "_"

			Plot_σΨm = FileSoilHydro_Plot * "/LabSigmaHm/" 
				mkpath(Plot_σΨm)
				path.plotSoilwater.Plot_σΨm = Plot_σΨm * SiteName * "_"

			Plot_Psd_θΨ = FileSoilHydro_Plot * "/Psd/IMP_ThetaH/"
				mkpath(Plot_Psd_θΨ)				
				path.plotSoilwater.Plot_Psd_θΨ = Plot_Psd_θΨ * SiteName * "_"

			Plot_IMP_model = FileSoilHydro_Plot * "/Psd/IMP/"
				mkpath(Plot_IMP_model)
				path.plotSoilwater.Plot_IMP_model = Plot_IMP_model * SiteName * "_"

			Plot_Psd_θr = FileSoilHydro_Plot * "/Psd/ThetaR/" 
				mkpath(Plot_Psd_θr)
				path.plotSoilwater.Plot_Psd_θr = Plot_Psd_θr * "Plot_ThetaR.svg"

			Plot_∑infilt_Opt = FileSoilHydro_Plot * "/Infiltration/Optimize/"
				mkpath(Plot_∑infilt_Opt)
				path.plotSoilwater.Plot_∑infilt_Opt = Plot_∑infilt_Opt * SiteName * "_"

			Plot_∑infilt_θΨ = FileSoilHydro_Plot * "/Infiltration/ThetaH/"
				mkpath(Plot_∑infilt_θΨ)
				path.plotSoilwater.Plot_∑infilt_θΨ = Plot_∑infilt_θΨ * SiteName * "_"
		
		# =============================================================
		#		HYPIX MODEL
		# =============================================================
			IdName_Hypix = Soilname[iSim]
		
		# HYPIX INPUT LEVEL 1 ===
         FileHypix_Input₁         = Path_Home * "/INPUT/Data_Hypix/" * SiteName * "/"
         path.hyPix.IdSelect      = FileHypix_Input₁ * SiteName * "_" * path.hyPix.IdSelect
         path.hyPix.JulesMetadata = FileHypix_Input₁ * path.hyPix.JulesMetadata

			# HYPIX INPUT LEVEL 2 ===
			FileHypix_Input₂            = Path_Home * "/INPUT/Data_Hypix/" * SiteName * "/" * IdName_Hypix * "/" * IdName_Hypix * "_"

			path.hyPix.Climate          = FileHypix_Input₂ * string(opt.hyPix.ClimateDataTimestep⍰) * "_" * path.hyPix.Climate
			path.hyPix.Dates            = FileHypix_Input₂ * path.hyPix.Dates
			path.hyPix.Discretization   = FileHypix_Input₂ * path.hyPix.Discretization
			path.hyPix.HyPix_HydroParam = FileHypix_Input₂ * path.hyPix.HyPix_HydroParam
			path.hyPix.HyPix_VegParam   = FileHypix_Input₂ * path.hyPix.HyPix_VegParam
			path.hyPix.HyPixParamOpt    = FileHypix_Input₂ * path.hyPix.HyPixParamOpt
			path.hyPix.obsTheta         = FileHypix_Input₂ * path.hyPix.obsTheta

			path.hyPix.Input_OfStep     = Path_Home * "/INPUT/Data_Hypix/RESULTS/"

			# HYPIX LOOKUPTABLE ===
         FileHypix_LookUpTable                 = Path_Home * "/INPUT/Data_Hypix/LookUpTable/"
					
         path.hyPix.LookUpTable_CropCoeficient = FileHypix_LookUpTable * path.hyPix.LookUpTable_CropCoeficient
         path.hyPix.LookUpTable_Lai                       = FileHypix_LookUpTable * path.hyPix.LookUpTable_Lai

			# HYPIX OUTPUT TABLE
			FileSoilHydro_Table = Path_Home * "/OUTPUT/Hypix/" * SiteName * "/" * IdName_Hypix *"/Table/" 				
				mkpath(FileSoilHydro_Table) #Make Folder if not exist

				FileSoilHydro_Table = FileSoilHydro_Table * SiteName * "_"

            path.hyPix.Table_DailyClimate    = FileSoilHydro_Table  *  IdName_Hypix * "_"* path.hyPix.Table_DailyClimate
            path.hyPix.Table_Discretisation  = FileSoilHydro_Table  *  IdName_Hypix * "_" *path.hyPix.Table_Discretisation
            path.hyPix.Table_Hydro           = FileSoilHydro_Table  *  IdName_Hypix * "_" *path.hyPix.Table_Hydro
            path.hyPix.Table_KΨ              = FileSoilHydro_Table  *  IdName_Hypix * "_"* path.hyPix.Table_KΨ
            path.hyPix.Table_Performance     = FileSoilHydro_Table  *  IdName_Hypix * "_"* path.hyPix.Table_Performance
            path.hyPix.Table_Q               = FileSoilHydro_Table  *  IdName_Hypix * "_"* path.hyPix.Table_Q
            path.hyPix.Table_Signature       = FileSoilHydro_Table  *  IdName_Hypix * "_"* path.hyPix.Table_Signature
            path.hyPix.Table_TimeSerie       = FileSoilHydro_Table  *  IdName_Hypix * "_"* path.hyPix.Table_TimeSerie
            path.hyPix.Table_TimeSerie_Daily = FileSoilHydro_Table  *  IdName_Hypix * "_"* path.hyPix.Table_TimeSerie_Daily
            path.hyPix.Table_Veg             = FileSoilHydro_Table  *  IdName_Hypix * "_"* path.hyPix.Table_Veg
            path.hyPix.Table_θ               = FileSoilHydro_Table  *  IdName_Hypix * "_"* path.hyPix.Table_θ
            path.hyPix.Table_θΨ              = FileSoilHydro_Table  *  IdName_Hypix * "_"* path.hyPix.Table_θΨ
            path.hyPix.Table_Ψ               = FileSoilHydro_Table  *  IdName_Hypix * "_"* path.hyPix.Table_Ψ
				
			FileSoilHydro_Table_θaverage = Path_Home * "/OUTPUT/Hypix/" * SiteName * "/SoilMoistureSim/" 				
				mkpath(FileSoilHydro_Table_θaverage) #Make Folder if not exist
				path.hyPix.Table_θaverage        = FileSoilHydro_Table_θaverage *  IdName_Hypix * "_"* path.hyPix.Table_θaverage

			# HYPIX PLOT CORE
			FileHypix_Plot = Path_Home * "/OUTPUT/Hypix/" * SiteName * "/" * IdName_Hypix *"/Plots/" 	
				mkpath(FileHypix_Plot)

            path.hyPix.Plot_HypixTime  = FileHypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_HypixTime
            path.hyPix.Plot_Hypix_θΨK  = FileHypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Hypix_θΨK
            path.hyPix.Plot_RainfallInterception  = FileHypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_RainfallInterception
            path.hyPix.Plot_Se_Time    = FileHypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Se_Time
            path.hyPix.Plot_Se_Z       = FileHypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Se_Z
            path.hyPix.Plot_Sorptivity = FileHypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Sorptivity
            path.hyPix.Vegetation      = FileHypix_Plot * IdName_Hypix  * "_" * path.hyPix.Vegetation

			# HYPIX PLOT OTHERS: RESULTS
			FileHypix_Plot_Results = Path_Home * "/OUTPUT/Hypix/RESULTS/"
				mkpath(FileHypix_Plot_Results)

		path.hyPix.Plot_OfStep           = FileHypix_Plot_Results
		path.hyPix.Plot_θ∂θ∂Ψ            = FileHypix_Plot_Results * path.hyPix.Plot_θ∂θ∂Ψ
		path.hyPix.Plot_Ψmin_Ψmax        = FileHypix_Plot_Results * path.hyPix.Plot_Ψmin_Ψmax
		path.hyPix.Plot_σ2θr             = FileHypix_Plot_Results * path.hyPix.Plot_σ2θr
		path.hyPix.Plot_θΨ_Δθ            = FileHypix_Plot_Results * path.hyPix.Plot_θΨ_Δθ
		path.hyPix.Plot_Se_Ψ_Constrained = FileHypix_Plot_Results * path.hyPix.Plot_Se_Ψ_Constrained

	
	return path
	end # function PATHS			
end  # module path
# ............................................................