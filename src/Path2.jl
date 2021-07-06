# =============================================================
#		MODULE: path
# =============================================================
module paths
	import ..tool
	import TOML

	mutable struct INPUT_SMAP
		LookupTable_RockWetability::String
		Smap::String
	end

	mutable struct INPUT_SOILWATER
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

	mutable struct INPUT_TEMPORARY
		σ_ψM_Scenario::String
	end
	
	struct OPTIONS 
		ModelName::String
		Select::String
	end # struct OPTION
	
	mutable struct PLOT_SOILWATER	
		Plot_∑infilt_Opt::String
		Plot_∑infilt_θΨ::String
		Plot_IMP_model::String
		Plot_Psd_θr::String
		Plot_Psd_θΨ::String				
		Plot_Sorptivity_Se::String
		Plot_θΨK::String
		Plot_σΨm::String
	end # PLOT_SOILWATER

	mutable struct TABLE_SMAP
		Table_Smap::String
		Table_θΨK::String
	end # struct INPUT_TABLE_SMAP

	mutable struct TABLE_SOILWATER
		FileSoilHydro_Table₁::String
		Table_HydroInfilt::String
		Table_Infilt::String
		Table_KΨ::String
		Table_Psd::String
		Table_Psd_θΨ_θ::String
		Table_θΨ_Psd::String
		Table_θΨK::String
		TableComplete_θΨ::String
		TableComplete_KΨ::String
	end # struct INPUT_SOILWATER

	mutable struct SMAP_2_HYPIX
		Path_Smap2Hypix::String
	end

	# _______________________ START: hypix _______________________ 
	mutable struct PATHYPIXS
		Climate::String
		Dates::String
		Discretization::String
		HyPix_HydroParam::String
		HyPixParamOpt::String

		HyPix_VegParam::String
		IdSelect::String 
		Input_OfStep::String
		JulesMetadata::String
		SiteName_Hypix::String
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

	mutable struct PATHS
		PathHome
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

		# _______________________ START: optionMaster _______________________
		struct OPTION_MASTER
			SiteName_Soilwater::String
		end 


	


		# ------------------------END: optionMaster---------------------------  
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PATHS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PATH(iSim, opt, optionMaster; Soilname=["Temporary"])
		PathHome₀ = @__DIR__

		# CHANGE DIRECTORY 
			PathHome₀ = dirname(PathHome₀)
			PathHome = PathHome₀ * "/data/"

		# READING OPTION_MASTER 
			Path_OptionMaster = PathHome₀ * "/data/OptionMaster.toml"
			Toml_OptionMaster = TOML.parsefile(Path_OptionMaster)

			optionMaster = OPTION_MASTER("SiteName_Soilwater")
			optionMaster = tool.readWrite.TOML_2_STRUCT(optionMaster, Toml_OptionMaster)

		

			
		# =============================================================
		#		OPTIONS
		# =============================================================
			ModelName ="Check"
			Select = "SELECT_1" # "SELECT_1" ;"SELECT_2" Select data to model

			option = OPTIONS(ModelName, Select)
	
			FileDataSoilhydro_Input = PathHome * "INPUT/Data_SoilWater/" * optionMaster.SiteName_Soilwater * "/" * optionMaster.SiteName_Soilwater * "_"

			Path =FileDataSoilhydro_Input * "Path.toml"
	
	
		# READING INPUT PATH.TOML
			TomlParse = TOML.parsefile(Path)


		# Paths
			FileDataSoilhydro_Input = PathHome * "INPUT/Data_SoilWater/" * optionMaster.SiteName_Soilwater * "/" * optionMaster.SiteName_Soilwater * "_"

			FileSoilHydro_Table₁ = PathHome * "/OUTPUT/SoilWater/" * optionMaster.SiteName_Soilwater * "/Table/" 
				mkpath(FileSoilHydro_Table₁) 

		# =============================================================
		#		INPUT_SOILWATER
		# =============================================================
			# DATA input into SoilWater
            # BulkDensity        = "BulkDensity.csv"
            # ConvertModel       = "TableHydro_Compiled_Homogeneous.csv"
            # HydroParam_Infilt  = "HydroParam_Infilt.csv"
            # HydroParam_ThetaH  = "GUI_HydroParam.csv"
            # IdSelect           = "IdSelect.csv"
            # Infiltration       = "Infiltration.csv"
            # Infiltration_Param = "InfiltrationParam.csv"
            # Kunsat             = "KunsatH.csv"
            # Psd                = "Psd.csv"
            # Pedological⍰      = "Pedological.csv"
            # Φ                  = "TotalPorosity.csv"
            # Ψθ                 = "ThetaH.csv"

				inputSoilwater = INPUT_SOILWATER("BulkDensity", "ConvertModel", "HydroParam_Infilt", "HydroParam_ThetaH", "IdSelect", "Infiltration", "Infiltration_Param", "Kunsat", "Psd", "Pedological⍰", "Φ", "Ψθ")
				
				inputSoilwater = tool.readWrite.TOML_2_STRUCT(inputSoilwater, TomlParse)	
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			inputSoilwater.BulkDensity        = FileDataSoilhydro_Input * inputSoilwater.inputSoilwater.BulkDensity
			inputSoilwater.ConvertModel       = FileDataSoilhydro_Input * inputSoilwater.ConvertModel
			inputSoilwater.HydroParam_Infilt  = FileDataSoilhydro_Input * inputSoilwater.HydroParam_Infilt
			inputSoilwater.HydroParam_ThetaH  = FileDataSoilhydro_Input * inputSoilwater.HydroParam_ThetaH
			inputSoilwater.IdSelect           = FileDataSoilhydro_Input * inputSoilwater.IdSelect
			inputSoilwater.Infiltration       = FileDataSoilhydro_Input * inputSoilwater.Infiltration
			inputSoilwater.Infiltration_Param = FileDataSoilhydro_Input * inputSoilwater.Infiltration_Param
			inputSoilwater.Kunsat             = FileDataSoilhydro_Input * inputSoilwater.Kunsat
			inputSoilwater.Psd                = FileDataSoilhydro_Input * inputSoilwater.Psd
			inputSoilwater.Pedological⍰    	= FileDataSoilhydro_Input * inputSoilwater.Pedological⍰
			inputSoilwater.Φ                  = FileDataSoilhydro_Input * inputSoilwater.Φ
			inputSoilwater.Ψθ                 = FileDataSoilhydro_Input * inputSoilwater.Ψθ

	
		# =============================================================
		#		INPUT_SMAP
		#		None core
		# =============================================================
			# Smap input path
				# LookupTable_RockWetability = "LookupTable_RockWetability.csv"				
				# Smap                    	= "SmapLayer.csv"
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			inputSmap = INPUT_SMAP("LookupTable_RockWetability", "Smap")
			inputSmap = tool.readWrite.TOML_2_STRUCT(inputSmap, TomlParse)	

         inputSmap.LookupTable_RockWetability = FileDataSoilhydro_Input * inputSmap.LookupTable_RockWetability
         inputSmap.Smap                       = FileDataSoilhydro_Input * inputSmap.Smap

			
				

		# =============================================================
		#		INPUT_TEMPORARY
		#		None core
		# =============================================================
			# Temporary input files	
				# σ_ψM_Scenario           = "σ_ψM_Scenario.csv"
			
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			inputTemporary = INPUT_TEMPORARY("σ_ψM_Scenario")
			inputTemporary = tool.readWrite.TOML_2_STRUCT(inputTemporary, TomlParse)	

			inputTemporary.σ_ψM_Scenario      = FileDataSoilhydro_Input * inputTemporary.σ_ψM_Scenario
				
		
		# =============================================================
		#		TABLE_SOILWATER
		# =============================================================
			# Output tables soil water
            # Table_KΨ            = "Table_KΨ.csv"
            # Table_HydroInfilt   = "Table_HydroInfilt.csv"
            # Table_Infilt        = "Table_Infilt.csv"
            # Table_Psd           = "Table_Psd.csv"
            # Table_Psd_θΨ_θ      = "Table_PsdTheta.csv"
            # Table_θΨ_Psd        = "Table_PsdHydro.csv"
            # Table_θΨK           = "Table_θΨK.csv"
            # TableComplete_θΨ    = "TableComplete_θΨ.csv"
				# TableComplete_KΨ    = "TableCompleteKΨ.csv"

				tableSoilwater = TABLE_SOILWATER("FileSoilHydro_Table₁", "Table_HydroInfilt", "Table_Infilt", "Table_KΨ", "Table_Psd", "Table_Psd_θΨ_θ", "Table_θΨ_Psd", "Table_θΨK", "TableComplete_θΨ", "TableComplete_KΨ")
				tableSoilwater = tool.readWrite.TOML_2_STRUCT(tableSoilwater, TomlParse)	
			
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			tableSoilwater.FileSoilHydro_Table₁ = tableSoilwater.FileSoilHydro_Table₁ * optionMaster.SiteName_Soilwater
			FileSoilHydro_Table₁ = FileSoilHydro_Table₁ * optionMaster.SiteName_Soilwater  

            tableSoilwater.Table_HydroInfilt    = FileSoilHydro_Table₁ * string(opt.infilt.Model⍰) * "_" *  ModelName  *  "_" * tableSoilwater.Table_HydroInfilt
            tableSoilwater.Table_Infilt         = FileSoilHydro_Table₁ * string(opt.infilt.Model⍰) *  "_" *  ModelName  *  "_" *  tableSoilwater.Table_Infilt
            tableSoilwater.Table_KΨ             = FileSoilHydro_Table₁ * "_"  *  ModelName * "_" * tableSoilwater.Table_KΨ
            tableSoilwater.Table_Psd            = FileSoilHydro_Table₁ * string(opt.psd.Model⍰) *  "_" * ModelName * "_" * tableSoilwater.Table_Psd
            tableSoilwater.Table_Psd_θΨ_θ       = FileSoilHydro_Table₁ * string(opt.psd.HydroModel⍰) *  "_" * ModelName * "_" *  tableSoilwater.Table_Psd_θΨ_θ
            tableSoilwater.Table_θΨ_Psd         = FileSoilHydro_Table₁ * string(opt.psd.HydroModel⍰) *  "_" * string(opt.hydro.σ_2_Ψm⍰) *  "_" * ModelName * "_" * tableSoilwater.Table_θΨ_Psd
            tableSoilwater.Table_θΨK            = FileSoilHydro_Table₁ *  "_" * string(opt.hydro.HydroModel⍰) * "_" * tableSoilwater.Table_θΨK
            tableSoilwater.TableComplete_θΨ     = FileSoilHydro_Table₁ *   "_" *  ModelName * "_" * tableSoilwater.TableComplete_θΨ
            tableSoilwater.TableComplete_KΨ     = FileSoilHydro_Table₁ *   "_" *  ModelName * "_" * tableSoilwater.TableComplete_KΨ
			
	
		# =============================================================
		#		TABLE_SMAP
		# =============================================================
			# Output tables smap
            # Table_θΨK = "Table_SmapThetaHK.csv"

				tableSmap = TABLE_SMAP("Table_Smap", "Table_θΨK")
				tableSmap = tool.readWrite.TOML_2_STRUCT(tableSmap, TomlParse)	

				tableSmap.Table_Smap = FileSoilHydro_Table₁ * "Table_Smap.csv"
            tableSmap.Table_θΨK = FileSoilHydro_Table₁ * "_" * string(opt.hydro.HydroModel⍰) * "_" * tableSmap.tableSmapTable_θΨK

				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		# =============================================================
		#		PATH SMAP_2_HYPIX
		# =============================================================
			# Output tables smap

			Path_Smap2Hypix = PathHome *"OUTPUT/Smap2Hypix"
			mkpath(Path_Smap2Hypix)

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		smap2Hypix = SMAP_2_HYPIX(Path_Smap2Hypix)

	
		# =============================================================
		#		PLOT SOILWATER
		# =============================================================
			FileSoilHydro_Plot = PathHome * "/OUTPUT/SoilWater/" * optionMaster.SiteName_Soilwater * "/Plots/"

			Plot_θΨK = FileSoilHydro_Plot * "/Lab/" 
				mkpath(Plot_θΨK)
				Plot_θΨK = Plot_θΨK * optionMaster.SiteName_Soilwater * "_"

			Plot_σΨm = FileSoilHydro_Plot * "/LabSigmaHm/" 
				mkpath(Plot_σΨm)
				Plot_σΨm = Plot_σΨm * optionMaster.SiteName_Soilwater * "_"

			Plot_Psd_θΨ = FileSoilHydro_Plot * "/Psd/IMP_ThetaH/"
				mkpath(Plot_Psd_θΨ)				
				Plot_Psd_θΨ = Plot_Psd_θΨ * optionMaster.SiteName_Soilwater * "_"

			Plot_IMP_model = FileSoilHydro_Plot * "/Psd/IMP/"
				mkpath(Plot_IMP_model)
				Plot_IMP_model = Plot_IMP_model * optionMaster.SiteName_Soilwater * "_"

			Plot_Psd_θr = FileSoilHydro_Plot * "/Psd/ThetaR/" 
				mkpath(Plot_Psd_θr)
				Plot_Psd_θr = Plot_Psd_θr * "Plot_ThetaR.svg"

			Plot_∑infilt_Opt = FileSoilHydro_Plot * "/Infiltration/Optimize/"
				mkpath(Plot_∑infilt_Opt)
				Plot_∑infilt_Opt = Plot_∑infilt_Opt * optionMaster.SiteName_Soilwater * "_"

			Plot_∑infilt_θΨ = FileSoilHydro_Plot * "/Infiltration/ThetaH/"
				mkpath(Plot_∑infilt_θΨ)
				Plot_∑infilt_θΨ = Plot_∑infilt_θΨ * optionMaster.SiteName_Soilwater * "_"

			Plot_Sorptivity_Se = FileSoilHydro_Plot * "/Sorptivity/"
				mkpath(Plot_Sorptivity_Se)
				Plot_Sorptivity_Se = Plot_Sorptivity_Se * optionMaster.SiteName_Soilwater * "_"

		plotSoilwater = PLOT_SOILWATER(Plot_∑infilt_Opt, Plot_∑infilt_θΨ, Plot_IMP_model, Plot_Psd_θr, Plot_Psd_θΨ, Plot_Sorptivity_Se, Plot_θΨK, Plot_σΨm)
		
		# =============================================================
		#		HYPIX MODEL
		# =============================================================
			# INPUT NAME OF FILE
				SiteName_Hypix = "LYSIMETERS" # "JULES"; "LYSIMETERS" 
			
				# TAUPO"; "OTOROHANGA"; "WAIHOU"; "WAITOA"; "HAMILTON"; 

		
			# # HYPIX INPUT JULES	
			# 	JulesMetadata = "JULES_LinkingData.csv"

			# # HYPIX INPUT DATA
         #    Climate          = "Climate_2.csv"
         #    Dates            = "Dates.csv"
         #    Discretization   = "Discretization_2.csv"
         #    HyPix_HydroParam = "HypixHydro.csv"
         #    HyPixParamOpt    = "HyPixParamOpt.csv"
         #    HyPix_VegParam   = "Vegetation.csv"
         #    IdSelect         = "IdSelect.csv"
         #    Input_OfStep     = "Wof_Steps.csv"
         #    obsTheta             = "Soilmoisture.csv"
				
			# # HYPIX LOOKUPTABLE
			# 	LookUpTable_CropCoeficient = "LookUpTable_CropCoeficient.csv"
			# 	LookUpTable_Lai            = "LookUpTable_Lai.csv"

			# # HYPIX OUTPUT TABLE
			# 	Table_DailyClimate    = "Table_DailyClimate"
			# 	Table_Discretisation  = "Table_Discretisation.csv"
			# 	Table_Hydro           = "Table_Hydro"
			# 	Table_KΨ              = "Table_KΨ"
			# 	Table_Performance     = "Table_Performance"
			# 	Table_Q               = "Table_Q"
			# 	Table_Signature       = "Table_Signature"
			# 	Table_TimeSerie       = "Table_TimeSerie"
			# 	Table_TimeSerie_Daily = "Table_TimeSerie_Daily"
			# 	Table_Veg             = "Table_Veg"
			# 	Table_Ψ               = "Table_H"
			# 	Table_θ               = "Table_Sm"
			# 	Table_θaverage        = "Table_THETAaverage"
			# 	Table_θΨ              = "Table_θΨ"

			# # HYPIX PLOTS 
			# 	Plot_HypixTime            = "Plot_HypixTime"
			# 	Plot_Hypix_θΨK            = "Plot_ThetaPsiK"
			# 	Plot_RainfallInterception = "Plot_RainfallInterception"
			# 	Plot_Se_Time              = "Plot_Se_Time.png"
			# 	Plot_Se_Z                 = "Plot_Se_Z.png"
			# 	Plot_Sorptivity           = "Plot_Sorptivity"
			# 	Vegetation                = "Plot_Vegetation"

			# # HYPIX PLOT OTHERS: RESULTS
			# 	# Plot_OfStep
         #    Plot_Se_Ψ_Constrained = "Plot_Se_Ψ_Constrained.svg"
         #    Plot_Ψmin_Ψmax        = "Plot_ΨminΨmax.svg"
         #    Plot_θΨ_Δθ            = "Plot_θΨ_Δθ.svg"
         #    Plot_σ2θr             = "Plot_θr2σ.svg"
         #    Plot_θ∂θ∂Ψ            = "Plot_θ∂θ∂Ψ.svg"
			
	
		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
		# 						PROCESSING DATA
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			IdName_Hypix = Soilname[iSim]

			# HYPIX INPUT JULES	
				JulesMetadata = "JULES_LinkingData.csv"

				hyPix = PATHYPIXS(
					"Climate",
					"Dates",
					"Discretization",
					"HyPix_HydroParam",
					"HyPixParamOpt",
					"HyPix_VegParam",
					"IdSelect",
					"Input_OfStep",
					"JulesMetadata",
					"SiteName_Hypix",
					"IdName_Hypix",
					"obsTheta", 
					"LookUpTable_CropCoeficient",
					"LookUpTable_Lai",
					"Table_DailyClimate",
					"Table_Discretisation",
					"Table_Hydro",
					"Table_KΨ",
					"Table_Performance",
					"Table_Q",
					"Table_Signature",
					"Table_TimeSerie",
					"Table_TimeSerie_Daily",
					"Table_Veg",
					"Table_Ψ",
					"Table_θ",
					"Table_θaverage",
					"Table_θΨ",
					"Plot_Hypix_θΨK",
					"Plot_HypixTime",
					"Plot_RainfallInterception",
					"Plot_Se_Time",
					"Plot_Se_Z",
					"Plot_Sorptivity",
					"Vegetation",
					"Plot_OfStep",
					"Plot_Se_Ψ_Constrained",
					"Plot_θΨ_Δθ",
					"Plot_σ2θr",
					"Plot_Ψmin_Ψmax",
					"Plot_θ∂θ∂Ψ")

				hyPix = tool.readWrite.TOML_2_STRUCT(hyPix, TomlParse)	

			# HYPIX INPUT LEVEL 1 ===
				FileHypix_Input₁ = PathHome * "/INPUT/Data_Hypix/" * SiteName_Hypix * "/"
				IdSelect =  FileHypix_Input₁ * SiteName_Hypix * "_" * IdSelect
				JulesMetadata   =FileHypix_Input₁ * JulesMetadata

			# HYPIX INPUT LEVEL 2 ===
				FileHypix_Input₂  = PathHome * "/INPUT/Data_Hypix/" * SiteName_Hypix * "/" * IdName_Hypix * "/" * IdName_Hypix * "_"

				hyPix.Climate          = FileHypix_Input₂ * string(opt.hyPix.ClimateDataTimestep⍰) * "_" * hyPix.Climate
				hyPix.Dates            = FileHypix_Input₂ * hyPix.Dates
				hyPix.Discretization   = FileHypix_Input₂ * hyPix.Discretization
				hyPix.HyPix_HydroParam = FileHypix_Input₂ * hyPix.HyPix_HydroParam
				hyPix.HyPix_VegParam   = FileHypix_Input₂ * hyPix.HyPix_VegParam
				hyPix.HyPixParamOpt      = FileHypix_Input₂ * hyPix.HyPixParamOpt
				hyPix.obsTheta             = FileHypix_Input₂ * hyPix.obsTheta

				hyPix.Input_OfStep     = PathHome * "/INPUT/Data_Hypix/RESULTS/"

			# HYPIX LOOKUPTABLE ===
			FileHypix_LookUpTable = PathHome * "/INPUT/Data_Hypix/LookUpTable/"
					
			hyPix.LookUpTable_CropCoeficient = FileHypix_LookUpTable * hyPix.LookUpTable_CropCoeficient
			hyPix.LookUpTable_Lai            = FileHypix_LookUpTable * hyPix.LookUpTable_Lai

			# HYPIX OUTPUT TABLE
			FileSoilHydro_Table = PathHome * "/OUTPUT/Hypix/" * SiteName_Hypix * "/" * IdName_Hypix *"/Table/" 				
				mkpath(FileSoilHydro_Table) #Make Folder if not exist

				FileSoilHydro_Table = FileSoilHydro_Table * SiteName_Hypix * "_"

				hyPix.Table_DailyClimate    = FileSoilHydro_Table  *  IdName_Hypix * "_"* hyPix.Table_DailyClimate
				hyPix.Table_Discretisation  = FileSoilHydro_Table  *  IdName_Hypix * "_" *hyPix.Table_Discretisation
				hyPix.Table_Hydro           = FileSoilHydro_Table  *  IdName_Hypix * "_" *hyPix.Table_Hydro
				hyPix.Table_KΨ              = FileSoilHydro_Table  *  IdName_Hypix * "_"* hyPix.Table_KΨ
				hyPix.Table_Performance     = FileSoilHydro_Table  *  IdName_Hypix * "_"* hyPix.Table_Performance
				hyPix.Table_Q               = FileSoilHydro_Table  *  IdName_Hypix * "_"* hyPix.Table_Q
				hyPix.Table_Signature       = FileSoilHydro_Table  *  IdName_Hypix * "_"* hyPix.Table_Signature
				hyPix.Table_TimeSerie       = FileSoilHydro_Table  *  IdName_Hypix * "_"* hyPix.Table_TimeSerie
				hyPix.Table_TimeSerie_Daily = FileSoilHydro_Table  *  IdName_Hypix * "_"* hyPix.Table_TimeSerie_Daily
				hyPix.Table_Veg             = FileSoilHydro_Table  *  IdName_Hypix * "_"* hyPix.Table_Veg
				hyPix.Table_θ               = FileSoilHydro_Table  *  IdName_Hypix * "_"* hyPix.Table_θ
				hyPix.Table_θΨ              = FileSoilHydro_Table  *  IdName_Hypix * "_"* hyPix.Table_θΨ
				hyPix.Table_Ψ               = FileSoilHydro_Table  *  IdName_Hypix * "_"* hyPix.Table_Ψ
				
			FileSoilHydro_Table_θaverage = PathHome * "/OUTPUT/Hypix/" * SiteName_Hypix * "/SoilMoistureSim/" 				
				mkpath(FileSoilHydro_Table_θaverage) #Make Folder if not exist
				Table_θaverage        = FileSoilHydro_Table_θaverage *  IdName_Hypix * "_"* Table_θaverage

			# FileSoilHydro_Table_Performace = PathHome * "/OUTPUT/Hypix/" * SiteName_Hypix * "/"
			# Table_Performance     = FileSoilHydro_Table_Performace  *  IdName_Hypix * "_"* string(iSim) * "_"* Table_Performance		
				


			# HYPIX PLOT CORE
			# FileHypix_Plot = PathHome * "/OUTPUT/Hypix/" * SiteName_Hypix * "/Plots/" 
			FileHypix_Plot = PathHome * "/OUTPUT/Hypix/" * SiteName_Hypix * "/" * IdName_Hypix *"/Plots/" 	
				mkpath(FileHypix_Plot)

				hyPix.Plot_HypixTime            = FileHypix_Plot * IdName_Hypix  * "_" * hyPix.Plot_HypixTime
				hyPix.Plot_Hypix_θΨK            = FileHypix_Plot * IdName_Hypix  * "_" * hyPix.Plot_Hypix_θΨK
				hyPix.Plot_RainfallInterception = FileHypix_Plot * IdName_Hypix  * "_" * hyPix.Plot_RainfallInterception
				hyPix.Plot_Se_Time              = FileHypix_Plot * IdName_Hypix  * "_" * hyPix.Plot_Se_Time
				hyPix.Plot_Se_Z                 = FileHypix_Plot * IdName_Hypix  * "_" * hyPix.Plot_Se_Z
				hyPix.Plot_Sorptivity           = FileHypix_Plot * IdName_Hypix  * "_" * hyPix.Plot_Sorptivity
				hyPix.Vegetation                = FileHypix_Plot * IdName_Hypix  * "_" * hyPix.Vegetation

			# HYPIX PLOT OTHERS: RESULTS
			FileHypix_Plot_Results = PathHome * "/OUTPUT/Hypix/RESULTS/"
				mkpath(FileHypix_Plot_Results)

            hyPix.Plot_OfStep           = FileHypix_Plot_Results
            hyPix.Plot_θ∂θ∂Ψ            = FileHypix_Plot_Results * hyPix.Plot_θ∂θ∂Ψ
            hyPix.Plot_Ψmin_Ψmax        = FileHypix_Plot_Results * hyPix.Plot_Ψmin_Ψmax
            hyPix.Plot_σ2θr             = FileHypix_Plot_Results * hyPix.Plot_σ2θr
            hyPix.Plot_θΨ_Δθ            = FileHypix_Plot_Results * hyPix.Plot_θΨ_Δθ
            hyPix.Plot_Se_Ψ_Constrained = FileHypix_Plot_Results * hyPix.Plot_Se_Ψ_Constrained

			# STRUCTURE
		path = PATHS(PathHome, hyPix, inputSmap, inputSoilwater, inputTemporary, option, plotSoilwater, smap2Hypix, tableSmap, tableSoilwater)

	return path
	end # function PATHS			
end  # module path
# ............................................................