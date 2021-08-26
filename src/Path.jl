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
		IdSelect::String
		Infiltration::String
		Infiltration_Param::String
		Kunsat::String
		Psd::String
		Pedological⍰::String
		Φ::String
		Ψθ::String
	end # struct INPUT_SOILWATER

	@option mutable struct INPUT_GUISOILWATER
		GUI_HydroParam::String
		GUI_KsModel::String 
	end # struct INPUT_SOILWATER

	@option mutable struct INPUT_TEMPORARY
		σ_ψM_Scenario::String
	end

	@option mutable struct TABLE_SMAP
		Table_Smap::String
		Table_θΨK::String
	end # struct INPUT_TABLE_SMAP	

	@option mutable struct TABLE_SOILWATER
		Path_Soilwater_Table::String
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
		inputGuiSoilwater::INPUT_GUISOILWATER
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
	function PATH(iSim, opt, PathData_Hypix, PathData_SoilWater, SiteName_Hypix, SiteName_Soilwater, Soilwater_OR_Hypix⍰; Soilname=["Temporary"])

		if Soilwater_OR_Hypix⍰ == "SoilWater"
			PathToml = PathData_SoilWater * "/ParamOptionPath/" * SiteName_Soilwater * "_Path.toml"

		elseif Soilwater_OR_Hypix⍰ == "Hypix"
			PathToml = PathData_Hypix * "/ParamOptionPath/" * SiteName_Hypix * "_Path.toml"
		end

		# Reading toml -> path structure
			path = Configurations.from_toml(PATHS, PathToml)

		# Usefull path
			Path_Home₀ = @__DIR__

		# perform cs..
			Path_Home = dirname(Path_Home₀)
			Path_Home = Path_Home * "/data/"
			
		# =============================================================
		#		INPUT_SOILWATER
		# =============================================================
			Path_Soilwater_Data = Path_Home * "INPUT/Data_SoilWater/" * SiteName_Soilwater * "/" * SiteName_Soilwater * "_"
			
         path.inputSoilwater.BulkDensity        = Path_Soilwater_Data * path.inputSoilwater.BulkDensity
         path.inputSoilwater.ConvertModel       = Path_Soilwater_Data * path.inputSoilwater.ConvertModel
         path.inputSoilwater.HydroParam_Infilt  = Path_Soilwater_Data * path.inputSoilwater.HydroParam_Infilt
			
         path.inputSoilwater.IdSelect           = Path_Soilwater_Data * path.inputSoilwater.IdSelect
         path.inputSoilwater.Infiltration       = Path_Soilwater_Data * path.inputSoilwater.Infiltration
         path.inputSoilwater.Infiltration_Param = Path_Soilwater_Data * path.inputSoilwater.Infiltration_Param
         path.inputSoilwater.Kunsat             = Path_Soilwater_Data * path.inputSoilwater.Kunsat
         path.inputSoilwater.Psd                = Path_Soilwater_Data * path.inputSoilwater.Psd
         path.inputSoilwater.Pedological⍰       = Path_Soilwater_Data * path.inputSoilwater.Pedological⍰
         path.inputSoilwater.Φ                  = Path_Soilwater_Data * path.inputSoilwater.Φ
         path.inputSoilwater.Ψθ                 = Path_Soilwater_Data * path.inputSoilwater.Ψθ
			
		# =============================================================
		#		INPUT_GUISOILWATER
		# =============================================================
			path.inputGuiSoilwater.GUI_HydroParam  = PathData_SoilWater * "/ParamOptionPath/" * SiteName_Soilwater * "_" * path.inputGuiSoilwater.GUI_HydroParam
			path.inputGuiSoilwater.GUI_KsModel  = PathData_SoilWater * "/ParamOptionPath/" * SiteName_Soilwater * "_" * path.inputGuiSoilwater.GUI_KsModel
		
		
		# =============================================================
		#		INPUT_SMAP
		# =============================================================
         path.inputSmap.LookupTable_RockWetability = Path_Soilwater_Data * path.inputSmap.LookupTable_RockWetability
         path.inputSmap.Smap                       = Path_Soilwater_Data * path.inputSmap.Smap

		# =============================================================
		#		INPUT_TEMPORARY
		# =============================================================
         path.inputTemporary.σ_ψM_Scenario = Path_Soilwater_Data * path.inputTemporary.σ_ψM_Scenario

		# =============================================================
		#		TABLE_SOILWATER
		# =============================================================
		Path_Soilwater_Table    = Path_Home * "/OUTPUT/SoilWater/" * SiteName_Soilwater * "/Table/"
				mkpath(Path_Soilwater_Table) 
	
			Path_Soilwater_Table                     = Path_Soilwater_Table * SiteName_Soilwater

			path.tableSoilwater.Path_Soilwater_Table = Path_Soilwater_Table
			path.tableSoilwater.Table_HydroInfilt    = Path_Soilwater_Table * path.option.ModelName * "_" *string(opt.infilt.Model⍰) * "_" *  path.option.ModelName  *  "_" * path.tableSoilwater.Table_HydroInfilt
			path.tableSoilwater.Table_Infilt         = Path_Soilwater_Table * path.option.ModelName * "_" *string(opt.infilt.Model⍰) *  "_" *  path.option.ModelName  *  "_" *  path.tableSoilwater.Table_Infilt
			path.tableSoilwater.Table_KΨ             = Path_Soilwater_Table * "_"  *  path.option.ModelName * "_" * path.tableSoilwater.Table_KΨ
			path.tableSoilwater.Table_Psd            = Path_Soilwater_Table *  "_" * path.option.ModelName * "_" *string(opt.psd.Model⍰) *  "_" * path.option.ModelName * "_" * path.tableSoilwater.Table_Psd
			path.tableSoilwater.Table_Psd_θΨ_θ       = Path_Soilwater_Table *  "_" * path.option.ModelName * "_" *string(opt.psd.HydroModel⍰) *  "_" * path.option.ModelName * "_" *  path.tableSoilwater.Table_Psd_θΨ_θ
			path.tableSoilwater.Table_θΨ_Psd         = Path_Soilwater_Table *  "_" * path.option.ModelName * "_" *string(opt.psd.HydroModel⍰) *  "_" * string(opt.hydro.σ_2_Ψm⍰) *  "_" * path.option.ModelName * "_" * path.tableSoilwater.Table_θΨ_Psd
			path.tableSoilwater.Table_θΨK            = Path_Soilwater_Table *   "_" *  path.option.ModelName * "_" * string(opt.hydro.HydroModel⍰) * "_" * path.tableSoilwater.Table_θΨK
			path.tableSoilwater.TableComplete_θΨ     = Path_Soilwater_Table *   "_" *  path.option.ModelName * "_" * path.tableSoilwater.TableComplete_θΨ
			path.tableSoilwater.TableComplete_KΨ     = Path_Soilwater_Table *   "_" *  path.option.ModelName * "_" * path.tableSoilwater.TableComplete_KΨ

		# =============================================================
		#		TABLE_SMAP
		# =============================================================
         path.tableSmap.Table_θΨK  = Path_Soilwater_Table * "_" * string(opt.hydro.HydroModel⍰) * "_" *  path.tableSmap.Table_θΨK
         path.tableSmap.Table_Smap = Path_Soilwater_Table * "_" * path.tableSmap.Table_Smap

		# =============================================================
		#		PATH SMAP_2_HYPIX
		# =============================================================
			path.smap2Hypix.Path_Smap2Hypix = Path_Home *"OUTPUT/Smap2Hypix"
				mkpath(path.smap2Hypix.Path_Smap2Hypix)

		# =============================================================
		#		PLOT SOILWATER
		# =============================================================
			Path_Soilwater_Plot = Path_Home * "/OUTPUT/SoilWater/" * SiteName_Soilwater * "/Plots/"

			Plot_θΨK = Path_Soilwater_Plot * "/Lab/" 
				mkpath(Plot_θΨK)
				path.plotSoilwater.Plot_θΨK = Plot_θΨK * SiteName_Soilwater * "_"

			Plot_σΨm = Path_Soilwater_Plot * "/LabSigmaHm/" 
				mkpath(Plot_σΨm)
				path.plotSoilwater.Plot_σΨm = Plot_σΨm * SiteName_Soilwater * "_"

			Plot_Psd_θΨ = Path_Soilwater_Plot * "/Psd/IMP_ThetaH/"
				mkpath(Plot_Psd_θΨ)				
				path.plotSoilwater.Plot_Psd_θΨ = Plot_Psd_θΨ * SiteName_Soilwater * "_"

			Plot_IMP_model = Path_Soilwater_Plot * "/Psd/IMP/"
				mkpath(Plot_IMP_model)
				path.plotSoilwater.Plot_IMP_model = Plot_IMP_model * SiteName_Soilwater * "_"

			Plot_Psd_θr = Path_Soilwater_Plot * "/Psd/ThetaR/" 
				mkpath(Plot_Psd_θr)
				path.plotSoilwater.Plot_Psd_θr = Plot_Psd_θr * "Plot_ThetaR.svg"

			Plot_∑infilt_Opt = Path_Soilwater_Plot * "/Infiltration/Optimize/"
				mkpath(Plot_∑infilt_Opt)
				path.plotSoilwater.Plot_∑infilt_Opt = Plot_∑infilt_Opt * SiteName_Soilwater * "_"

			Plot_∑infilt_θΨ = Path_Soilwater_Plot * "/Infiltration/ThetaH/"
				mkpath(Plot_∑infilt_θΨ)
				path.plotSoilwater.Plot_∑infilt_θΨ = Plot_∑infilt_θΨ * SiteName_Soilwater * "_"
		
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

			path.hyPix.Climate          = Path_Hypix₂ * string(opt.hyPix.ClimateDataTimestep⍰) * "_" * path.hyPix.Climate
			path.hyPix.Dates            = Path_Hypix₂ * path.hyPix.Dates
			path.hyPix.Discretization   = Path_Hypix₂ * path.hyPix.Discretization
			path.hyPix.HyPix_HydroParam = Path_Hypix₂ * path.hyPix.HyPix_HydroParam
			path.hyPix.HyPix_VegParam   = Path_Hypix₂ * path.hyPix.HyPix_VegParam
			path.hyPix.HyPixParamOpt    = Path_Hypix₂ * path.hyPix.HyPixParamOpt
			path.hyPix.obsTheta         = Path_Hypix₂ * path.hyPix.obsTheta

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

            path.hyPix.Plot_HypixTime  = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_HypixTime
            path.hyPix.Plot_Hypix_θΨK  = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Hypix_θΨK
            path.hyPix.Plot_RainfallInterception  = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_RainfallInterception
            path.hyPix.Plot_Se_Time    = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Se_Time
            path.hyPix.Plot_Se_Z       = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Se_Z
            path.hyPix.Plot_Sorptivity = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Plot_Sorptivity
            path.hyPix.Vegetation      = Path_Hypix_Plot * IdName_Hypix  * "_" * path.hyPix.Vegetation

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