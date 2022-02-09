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
		HydroParamPrecomputed::String
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

	@option mutable struct  CONVERT_SOILWATER
   	Table_Convert_θΨ_2D_2_1D::String
   	Table_Convert_KΨ_2D_2_1D::String
	end

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
		Table_KsModel_Ks::String
		Table_KsModel::String
		Table_KsModel_τ::String
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
		Plot_KsModel::String
		Plot_Psd_θr::String
		Plot_Psd_θΨ::String				
		Plot_θΨK::String
		Plot_σΨm::String
	end # PLOT_SOILWATER

	@option mutable struct SMAP_2_HYPIX
		Path_Smap2Hypix::String
	end

	# ------------------------END: hypix---------------------------  
	@option mutable struct PATHS
		inputSmap::INPUT_SMAP
		inputSoilwater::INPUT_SOILWATER
		inputGuiSoilwater::INPUT_GUISOILWATER
		convertSoilwater::CONVERT_SOILWATER
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
			
         path.inputSoilwater.BulkDensity           = Path_Soilwater_Data * path.inputSoilwater.BulkDensity
         path.inputSoilwater.ConvertModel          = Path_Soilwater_Data * path.inputSoilwater.ConvertModel
         path.inputSoilwater.HydroParamPrecomputed = Path_Soilwater_Data * path.inputSoilwater.HydroParamPrecomputed
         path.inputSoilwater.HydroParam_Infilt     = Path_Soilwater_Data * path.inputSoilwater.HydroParam_Infilt

         path.inputSoilwater.IdSelect              = Path_Soilwater_Data * path.inputSoilwater.IdSelect
         path.inputSoilwater.Infiltration          = Path_Soilwater_Data * path.inputSoilwater.Infiltration
         path.inputSoilwater.Infiltration_Param    = Path_Soilwater_Data * path.inputSoilwater.Infiltration_Param
         path.inputSoilwater.Kunsat                = Path_Soilwater_Data * path.inputSoilwater.Kunsat
         path.inputSoilwater.Pedological⍰          = Path_Soilwater_Data * path.inputSoilwater.Pedological⍰
         path.inputSoilwater.Psd                   = Path_Soilwater_Data * path.inputSoilwater.Psd
         path.inputSoilwater.Φ                     = Path_Soilwater_Data * path.inputSoilwater.Φ
         path.inputSoilwater.Ψθ                    = Path_Soilwater_Data * path.inputSoilwater.Ψθ
			
		# =============================================================
		#		INPUT_GUISOILWATER
		# =============================================================
			path.inputGuiSoilwater.GUI_HydroParam  = PathData_SoilWater * "/ParamOptionPath/" * SiteName_Soilwater * "_" * path.inputGuiSoilwater.GUI_HydroParam
			path.inputGuiSoilwater.GUI_KsModel  = PathData_SoilWater * "/ParamOptionPath/" * SiteName_Soilwater * "_" * path.inputGuiSoilwater.GUI_KsModel

		# =============================================================
		#		CONVERT_SOILWATER
		# =============================================================
			Path_ConvertSoilWater = Path_Home * "INPUT/Data_SoilWater/" * SiteName_Soilwater * "/" * "Convert/"
			mkpath(Path_ConvertSoilWater)
				Path_ConvertSoilWater = Path_ConvertSoilWater * SiteName_Soilwater * "_"

				path.convertSoilwater.Table_Convert_θΨ_2D_2_1D = Path_ConvertSoilWater * path.convertSoilwater.Table_Convert_θΨ_2D_2_1D
				path.convertSoilwater.Table_Convert_KΨ_2D_2_1D = Path_ConvertSoilWater * path.convertSoilwater.Table_Convert_KΨ_2D_2_1D

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
			path.tableSoilwater.Table_KsModel_Ks        = Path_Soilwater_Table * "_"  *  path.option.ModelName * "_" *  path.option.Select  * "_" * path.tableSoilwater.Table_KsModel_Ks
			path.tableSoilwater.Table_KsModel        = Path_Soilwater_Table * "_"  *  path.option.ModelName * "_" *  opt.ksModel.KₛModel⍰ * "_" * path.tableSoilwater.Table_KsModel
			path.tableSoilwater.Table_KsModel_τ        = Path_Soilwater_Table * "_"  *  path.option.ModelName * "_" *  opt.ksModel.KₛModel⍰  * "_" * path.tableSoilwater.Table_KsModel_τ
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
				path.plotSoilwater.Plot_θΨK = Plot_θΨK * SiteName_Soilwater * "_" * path.option.ModelName *  "_"

			Plot_KsModel₀ = Path_Soilwater_Plot * "/KsModel/"
				mkpath(Plot_KsModel₀)
				path.plotSoilwater.Plot_KsModel = Plot_KsModel₀ * SiteName_Soilwater * "_"  * path.option.ModelName  * "_" * opt.ksModel.KₛModel⍰  * "_" * path.plotSoilwater.Plot_KsModel


			Plot_σΨm = Path_Soilwater_Plot * "/LabSigmaHm/" 
				mkpath(Plot_σΨm)
				path.plotSoilwater.Plot_σΨm = Plot_σΨm * SiteName_Soilwater * "_"  * path.option.ModelName *  "_"

			Plot_Psd_θΨ = Path_Soilwater_Plot * "/Psd/IMP_ThetaH/"
				mkpath(Plot_Psd_θΨ)				
				path.plotSoilwater.Plot_Psd_θΨ = Plot_Psd_θΨ * SiteName_Soilwater * "_"  * path.option.ModelName *  "_"

			Plot_IMP_model = Path_Soilwater_Plot * "/Psd/IMP/"
				mkpath(Plot_IMP_model)
				path.plotSoilwater.Plot_IMP_model = Plot_IMP_model * SiteName_Soilwater * "_"  * path.option.ModelName *  "_"

			Plot_Psd_θr = Path_Soilwater_Plot * "/Psd/ThetaR/" 
				mkpath(Plot_Psd_θr)
				path.plotSoilwater.Plot_Psd_θr = Plot_Psd_θr * "Plot_ThetaR.svg"

			Plot_∑infilt_Opt = Path_Soilwater_Plot * "/Infiltration/Optimize/"
				mkpath(Plot_∑infilt_Opt)
				path.plotSoilwater.Plot_∑infilt_Opt = Plot_∑infilt_Opt * SiteName_Soilwater * "_"  * path.option.ModelName *  "_"

			Plot_∑infilt_θΨ = Path_Soilwater_Plot * "/Infiltration/ThetaH/"
				mkpath(Plot_∑infilt_θΨ)
				path.plotSoilwater.Plot_∑infilt_θΨ = Plot_∑infilt_θΨ * SiteName_Soilwater * "_"  * path.option.ModelName *  "_"
	
	return path
	end # function PATHS			
end  # module path
# ............................................................