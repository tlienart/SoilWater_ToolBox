# =============================================================
#		MODULE: path
# =============================================================
module paths

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
		SiteName_Soilhyro::String
	end # struct OPTION
	
	mutable struct PLOT_SOILWATER	
		Plot_∑infilt_Opt::String
		Plot_∑infilt_SeIniRange::String
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
				HyPix_Param::String

				HyPix_VegParam::String
				Hydraulic_Kg::String
				IdSelect::String 
				Input_OfStep::String
				JulesMetadata::String
				ProjectName_Hypix::String
				IdName_Hypix::String
				obsθ::String 
		
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
			Home
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
	function PATH(iSim, opt; Soilname=["Temporary"])
		Home2 = @__DIR__

		# perform cs..
			Home = dirname(Home2)

		# Change path name to /data/Private/
			Home = Home * "/data/" * opt.other.DataPrivateShare * "/"
			
		# =============================================================
		#		OPTIONS
		# =============================================================
			# Which files to use
			SiteName_Soilhyro =  "NewFormat" #"Smap20210226"; "VCSNSmap2"; "SFF"; "PAF"; K10KPA; Smap; Smap20210226; SmapSouthland2; CantyLysimSmap; VCSNSmap; "WaikLysim"; "Convert; "SmapNZAllSoilsSmap20210326"; "Smap20210226"
			ModelName ="Check"
			Select = "SELECT_1" # Select data to model

			option = OPTIONS(ModelName, Select, SiteName_Soilhyro)
	
		# Paths
			FileDataSoilhydro_Input = Home * "INPUT/Data_SoilWater/" * SiteName_Soilhyro * "/" * SiteName_Soilhyro * "_"

			FileSoilHydro_Table₁ = Home * "/OUTPUT/SoilWater/" * SiteName_Soilhyro * "/Table/" 
				mkpath(FileSoilHydro_Table₁) 

		# =============================================================
		#		INPUT_SOILWATER
		# =============================================================
			# DATA input into SoilWater
            BulkDensity        = "BulkDensity.csv"
            ConvertModel       = "TableHydro_Compiled_Homogeneous.csv"
            HydroParam_Infilt  = "HydroParam_Infilt.csv"
            HydroParam_ThetaH  = "GUI_HydroParam.csv"
            IdSelect           = "IdSelect.csv"
            Infiltration       = "Infiltration.csv"
            Infiltration_Param = "InfiltrationParam.csv"
            Kunsat             = "KunsatH.csv"
            Psd                = "Psd.csv"
            Pedological⍰      = "Pedological.csv"
            Φ                  = "TotalPorosity.csv"
            Ψθ                 = "ThetaH.csv"
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            BulkDensity        = FileDataSoilhydro_Input * BulkDensity
            ConvertModel       = FileDataSoilhydro_Input * ConvertModel
            HydroParam_Infilt  = FileDataSoilhydro_Input * HydroParam_Infilt
            HydroParam_ThetaH  = FileDataSoilhydro_Input * HydroParam_ThetaH
            IdSelect           = FileDataSoilhydro_Input * IdSelect
            Infiltration       = FileDataSoilhydro_Input * Infiltration
            Infiltration_Param = FileDataSoilhydro_Input * Infiltration_Param
            Kunsat             = FileDataSoilhydro_Input * Kunsat
            Psd                = FileDataSoilhydro_Input * Psd
            Pedological⍰    	= FileDataSoilhydro_Input * Pedological⍰
            Φ                  = FileDataSoilhydro_Input * Φ
            Ψθ                 = FileDataSoilhydro_Input * Ψθ

			inputSoilwater = INPUT_SOILWATER(BulkDensity, ConvertModel, HydroParam_Infilt, HydroParam_ThetaH, IdSelect, Infiltration, Infiltration_Param, Kunsat, Psd, Pedological⍰, Φ, Ψθ)		
				
		# =============================================================
		#		INPUT_SMAP
		#		None core
		# =============================================================
			# Smap input path
				LookupTable_RockWetability = "LookupTable_RockWetability.csv"				
				Smap                    	= "SmapLayer.csv"
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            LookupTable_RockWetability = FileDataSoilhydro_Input * LookupTable_RockWetability
				Smap = FileDataSoilhydro_Input * Smap

			inputSmap = INPUT_SMAP(LookupTable_RockWetability, Smap)
				

		# =============================================================
		#		INPUT_TEMPORARY
		#		None core
		# =============================================================
			# Temporary input files	
				σ_ψM_Scenario           = "σ_ψM_Scenario.csv"
			
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				σ_ψM_Scenario      = FileDataSoilhydro_Input * σ_ψM_Scenario
				
			inputTemporary = INPUT_TEMPORARY(σ_ψM_Scenario)


		# =============================================================
		#		TABLE_SOILWATER
		# =============================================================
			# Output tables soil water
            Table_KΨ            = "Table_KΨ.csv"
            Table_HydroInfilt   = "Table_HydroInfilt.csv"
            Table_Infilt        = "Table_Infilt.csv"
            Table_Psd           = "Table_Psd.csv"
            Table_Psd_θΨ_θ      = "Table_PsdTheta.csv"
            Table_θΨ_Psd        = "Table_PsdHydro.csv"
            Table_θΨK           = "Table_θΨK.csv"
            TableComplete_θΨ    = "TableComplete_θΨ.csv"
				TableComplete_KΨ    = "TableCompleteKΨ.csv"
			
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            FileSoilHydro_Table₁ = FileSoilHydro_Table₁ * SiteName_Soilhyro

            Table_HydroInfilt    = FileSoilHydro_Table₁ * string(opt.infilt.Model) * "_" *  ModelName  *  "_" * Table_HydroInfilt
            Table_Infilt         = FileSoilHydro_Table₁ * string(opt.infilt.Model) *  "_" *  ModelName  *  "_" *  Table_Infilt
            Table_KΨ             = FileSoilHydro_Table₁ * "_"  *  ModelName * "_" * Table_KΨ
            Table_Psd            = FileSoilHydro_Table₁ * string(opt.psd.Model) *  "_" * ModelName * "_" * Table_Psd
            Table_Psd_θΨ_θ       = FileSoilHydro_Table₁ * string(opt.psd.HydroModel⍰) *  "_" * ModelName * "_" *  Table_Psd_θΨ_θ
            Table_θΨ_Psd         = FileSoilHydro_Table₁ * string(opt.psd.HydroModel⍰) *  "_" * string(opt.hydro.σ_2_Ψm⍰) *  "_" * ModelName * "_" * Table_θΨ_Psd
            Table_θΨK            = FileSoilHydro_Table₁ *  "_" * string(opt.hydro.HydroModel⍰) * "_" * Table_θΨK
            TableComplete_θΨ     = FileSoilHydro_Table₁ *   "_" *  ModelName * "_" * TableComplete_θΨ
            TableComplete_KΨ     = FileSoilHydro_Table₁ *   "_" *  ModelName * "_" * TableComplete_KΨ
			
			tableSoilwater = TABLE_SOILWATER(FileSoilHydro_Table₁, Table_HydroInfilt, Table_Infilt, Table_KΨ, Table_Psd, Table_Psd_θΨ_θ, Table_θΨ_Psd, Table_θΨK, TableComplete_θΨ, TableComplete_KΨ)


		# =============================================================
		#		TABLE_SMAP
		# =============================================================
			# Output tables smap
            Table_θΨK = "Table_SmapThetaHK.csv"

            Table_θΨK = FileSoilHydro_Table₁ * "_" * string(opt.hydro.HydroModel⍰) * "_" * Table_θΨK

				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				Table_Smap =  FileSoilHydro_Table₁ * "Table_Smap.csv"

				tableSmap = TABLE_SMAP(Table_Smap, Table_θΨK)

		# =============================================================
		#		PATH SMAP_2_HYPIX
		# =============================================================
			# Output tables smap

			Path_Smap2Hypix = Home *"OUTPUT/Smap2Hypix"
			mkpath(Path_Smap2Hypix)

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		smap2Hypix = SMAP_2_HYPIX(Path_Smap2Hypix)


		
		# =============================================================
		#		PLOT SOILWATER
		# =============================================================
			FileSoilHydro_Plot = Home * "/OUTPUT/SoilWater/" * SiteName_Soilhyro * "/Plots/"

			Plot_θΨK = FileSoilHydro_Plot * "/Lab/" 
				mkpath(Plot_θΨK)
				Plot_θΨK = Plot_θΨK * SiteName_Soilhyro * "_"

			Plot_σΨm = FileSoilHydro_Plot * "/LabSigmaHm/" 
				mkpath(Plot_σΨm)
				Plot_σΨm = Plot_σΨm * SiteName_Soilhyro * "_"

			Plot_Psd_θΨ = FileSoilHydro_Plot * "/Psd/IMP_ThetaH/"
				mkpath(Plot_Psd_θΨ)				
				Plot_Psd_θΨ = Plot_Psd_θΨ * SiteName_Soilhyro * "_"

			Plot_IMP_model = FileSoilHydro_Plot * "/Psd/IMP/"
				mkpath(Plot_IMP_model)
				Plot_IMP_model = Plot_IMP_model * SiteName_Soilhyro * "_"

			Plot_Psd_θr = FileSoilHydro_Plot * "/Psd/ThetaR/" 
				mkpath(Plot_Psd_θr)
				Plot_Psd_θr = Plot_Psd_θr * "Plot_ThetaR.svg"

			Plot_∑infilt_Opt = FileSoilHydro_Plot * "/Infiltration/Optimize/"
				mkpath(Plot_∑infilt_Opt)
				Plot_∑infilt_Opt = Plot_∑infilt_Opt * SiteName_Soilhyro * "_"

			Plot_∑infilt_SeIniRange = FileSoilHydro_Plot * "/Infiltration/SeIni/"
				mkpath(Plot_∑infilt_SeIniRange)
				Plot_∑infilt_SeIniRange = Plot_∑infilt_SeIniRange * SiteName_Soilhyro * "_"

			Plot_∑infilt_θΨ = FileSoilHydro_Plot * "/Infiltration/ThetaH/"
				mkpath(Plot_∑infilt_θΨ)
				Plot_∑infilt_θΨ = Plot_∑infilt_θΨ * SiteName_Soilhyro * "_"

			Plot_Sorptivity_Se = FileSoilHydro_Plot * "/Infiltration/Sorptivity/"
				mkpath(Plot_Sorptivity_Se)
				Plot_Sorptivity_Se = Plot_Sorptivity_Se * SiteName_Soilhyro * "_"

		plotSoilwater = PLOT_SOILWATER(Plot_∑infilt_Opt, Plot_∑infilt_SeIniRange, Plot_∑infilt_θΨ, Plot_IMP_model, Plot_Psd_θr, Plot_Psd_θΨ, Plot_Sorptivity_Se, Plot_θΨK, Plot_σΨm)
		
		# =============================================================
		#		HYPIX MODEL
		# =============================================================
			# INPUT NAME OF FILE
				ProjectName_Hypix = "JULES" # "JULES"; "LYSIMETERS" 
			
				# IdName_Hypix = "Lincoln" # "TAUPO"; "OTOROHANGA"; "WAIHOU"; "WAITOA"; "HAMILTON"; "Lincoln";
				IdName_Hypix = Soilname[iSim]
		
			# HYPIX INPUT JULES	
				JulesMetadata = "JULES_LinkingData.csv"

			# HYPIX INPUT DATA
            Climate          = "Climate_2.csv"
            Dates            = "Dates.csv"
            Discretization   = "Discretization_2.csv"
            HyPix_HydroParam = "HypixHydro.csv"
            HyPix_Param      = "HyPix_Param_2.csv"
            HyPix_VegParam   = "Vegetation.csv"
            Hydraulic_Kg     = "Hydraulic_Uni_Kg2.csv"
            IdSelect         = "IdSelect.csv"
            Input_OfStep     = "Wof_Steps.csv"
            obsθ             = "Soilmoisture.csv"
				
			# HYPIX LOOKUPTABLE
				LookUpTable_CropCoeficient = "LookUpTable_CropCoeficient.csv"
				LookUpTable_Lai            = "LookUpTable_Lai.csv"

			# HYPIX OUTPUT TABLE
				Table_DailyClimate    = "Table_DailyClimate"
				Table_Discretisation  = "Table_Discretisation.csv"
				Table_Hydro           = "Table_Hydro"
				Table_KΨ              = "Table_KΨ"
				Table_Performance     = "Table_Performance"
				Table_Q               = "Table_Q"
				Table_Signature       = "Table_Signature"
				Table_TimeSerie       = "Table_TimeSerie"
				Table_TimeSerie_Daily = "Table_TimeSerie_Daily"
				Table_Veg             = "Table_Veg"
				Table_Ψ               = "Table_H"
				Table_θ               = "Table_Sm"
				Table_θaverage        = "Table_THETAaverage"
				Table_θΨ              = "Table_θΨ"

			# HYPIX PLOTS 
				Plot_HypixTime            = "Plot_HypixTime"
				Plot_Hypix_θΨK            = "Plot_ThetaPsiK"
				Plot_RainfallInterception = "Plot_RainfallInterception"
				Plot_Se_Time              = "Plot_Se_Time.png"
				Plot_Se_Z                 = "Plot_Se_Z.png"
				Plot_Sorptivity           = "Plot_Sorptivity"
				Vegetation                = "Plot_Vegetation"

			# HYPIX PLOT OTHERS: RESULTS
				# Plot_OfStep
            Plot_Se_Ψ_Constrained = "Plot_Se_Ψ_Constrained.svg"
            Plot_Ψmin_Ψmax        = "Plot_ΨminΨmax.svg"
            Plot_θΨ_Δθ            = "Plot_θΨ_Δθ.svg"
            Plot_σ2θr             = "Plot_θr2σ.svg"
            Plot_θ∂θ∂Ψ            = "Plot_θ∂θ∂Ψ.svg"
	
		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
		# 						PROCESSING DATA
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			# INPUT NAME OF FILE
				ProjectName_Hypix = "JULES" # "JULES"; "LYSIMETERS" 
				
			# IdName_Hypix = "Lincoln" # "TAUPO"; "OTOROHANGA"; "WAIHOU"; "WAITOA"; "HAMILTON"; "Lincoln";
				IdName_Hypix = Soilname[iSim]
	
			# HYPIX INPUT JULES	
				JulesMetadata = "JULES_LinkingData.csv"

			# HYPIX INPUT LEVEL 1 ===
				FileHypix_Input₁ = Home * "/INPUT/Data_Hypix/" * ProjectName_Hypix * "/"
				IdSelect =  FileHypix_Input₁ * ProjectName_Hypix * "_" * IdSelect
				JulesMetadata   =FileHypix_Input₁ * JulesMetadata

			# HYPIX INPUT LEVEL 2 ===
				FileHypix_Input₂  = Home * "/INPUT/Data_Hypix/" * ProjectName_Hypix * "/" * IdName_Hypix * "/" * IdName_Hypix * "_"


				Climate          = FileHypix_Input₂ * opt.hyPix.ClimateDataTimestep * "_" * Climate
				Dates            = FileHypix_Input₂ * Dates
				Discretization   = FileHypix_Input₂ * Discretization
				HyPix_HydroParam = FileHypix_Input₂ * HyPix_HydroParam
				HyPix_VegParam   = FileHypix_Input₂ * HyPix_VegParam
				Hypix_Param      = FileHypix_Input₂ * HyPix_Param
				obsθ             = FileHypix_Input₂ * obsθ

				Input_OfStep     = Home * "/INPUT/Data_Hypix/RESULTS/"

			# HYPIX LOOKUPTABLE ===
			FileHypix_LookUpTable = Home * "/INPUT/Data_Hypix/LookUpTable/"
					
				LookUpTable_CropCoeficient = FileHypix_LookUpTable * LookUpTable_CropCoeficient
				LookUpTable_Lai            = FileHypix_LookUpTable * LookUpTable_Lai

			# HYPIX OUTPUT TABLE
			FileSoilHydro_Table = Home * "/OUTPUT/Hypix/" * ProjectName_Hypix * "/" * IdName_Hypix *"/Table/" 				
				mkpath(FileSoilHydro_Table) #Make Folder if not exist

				FileSoilHydro_Table = FileSoilHydro_Table * ProjectName_Hypix * "_"

				Table_DailyClimate    = FileSoilHydro_Table  *  IdName_Hypix * "_"* Table_DailyClimate
				Table_Discretisation  = FileSoilHydro_Table  *  IdName_Hypix * "_" *Table_Discretisation
				Table_Hydro           = FileSoilHydro_Table  *  IdName_Hypix * "_" *Table_Hydro
				Table_KΨ              = FileSoilHydro_Table  *  IdName_Hypix * "_"* Table_KΨ
				Table_Q               = FileSoilHydro_Table  *  IdName_Hypix * "_"* Table_Q
				Table_Signature       = FileSoilHydro_Table  *  IdName_Hypix * "_"* Table_Signature
				Table_TimeSerie       = FileSoilHydro_Table  *  IdName_Hypix * "_"* Table_TimeSerie
				Table_TimeSerie_Daily = FileSoilHydro_Table  *  IdName_Hypix * "_"* Table_TimeSerie_Daily
				Table_Veg             = FileSoilHydro_Table  *  IdName_Hypix * "_"* Table_Veg
				Table_Ψ               = FileSoilHydro_Table  *  IdName_Hypix * "_"* Table_Ψ
				Table_θ               = FileSoilHydro_Table  *  IdName_Hypix * "_"* Table_θ
				Table_θΨ              = FileSoilHydro_Table  *  IdName_Hypix * "_"* Table_θΨ
				
			FileSoilHydro_Table_θaverage = Home * "/OUTPUT/Hypix/" * ProjectName_Hypix * "/SoilMoistureSim/" 				
				mkpath(FileSoilHydro_Table_θaverage) #Make Folder if not exist
				Table_θaverage        = FileSoilHydro_Table_θaverage *  IdName_Hypix * "_"* Table_θaverage

			FileSoilHydro_Table_Performace = Home * "/OUTPUT/Hypix/" * ProjectName_Hypix * "/"		
				Table_Performance     = FileSoilHydro_Table_Performace  *  IdName_Hypix * "_"* string(iSim) * "_"* Table_Performance


			# HYPIX PLOT CORE
			FileHypix_Plot = Home * "/OUTPUT/Hypix/" * ProjectName_Hypix * "/Plots/" 
				mkpath(FileHypix_Plot)

				Plot_HypixTime            = FileHypix_Plot * IdName_Hypix  * "_" * Plot_HypixTime
				Plot_Hypix_θΨK            = FileHypix_Plot * IdName_Hypix  * "_" * Plot_Hypix_θΨK
				Plot_RainfallInterception = FileHypix_Plot * IdName_Hypix  * "_" * Plot_RainfallInterception
				Plot_Se_Time              = FileHypix_Plot * IdName_Hypix  * "_" * Plot_Se_Time
				Plot_Se_Z                 = FileHypix_Plot * IdName_Hypix  * "_" * Plot_Se_Z
				Plot_Sorptivity           = FileHypix_Plot * IdName_Hypix  * "_" * Plot_Sorptivity
				Vegetation                = FileHypix_Plot * IdName_Hypix  * "_" * Vegetation

			# HYPIX PLOT OTHERS: RESULTS
			FileHypix_Plot_Results = Home * "/OUTPUT/Hypix/RESULTS/"
				mkpath(FileHypix_Plot_Results)

				Plot_OfStep   = FileHypix_Plot_Results
				Plot_θ∂θ∂Ψ    = FileHypix_Plot_Results * Plot_θ∂θ∂Ψ
				Plot_Ψmin_Ψmax = FileHypix_Plot_Results * Plot_Ψmin_Ψmax
				Plot_σ2θr      = FileHypix_Plot_Results * Plot_σ2θr
				Plot_θΨ_Δθ     = FileHypix_Plot_Results * Plot_θΨ_Δθ
				Plot_Se_Ψ_Constrained = FileHypix_Plot_Results * Plot_Se_Ψ_Constrained

			# STRUCTURE
				hyPix = PATHYPIXS(
					Climate,
					Dates,
					Discretization,
					HyPix_HydroParam,
					HyPix_Param,
					HyPix_VegParam,
					Hydraulic_Kg,
					IdSelect,
					Input_OfStep,
					JulesMetadata,
					ProjectName_Hypix,
					IdName_Hypix,
					obsθ, 

					LookUpTable_CropCoeficient,
					LookUpTable_Lai,

					Table_DailyClimate,
					Table_Discretisation,
					Table_Hydro,
					Table_KΨ,
					Table_Performance,
					Table_Q,
					Table_Signature,
					Table_TimeSerie,
					Table_TimeSerie_Daily,
					Table_Veg,
					Table_Ψ,
					Table_θ,
					Table_θaverage,
					Table_θΨ,

					Plot_Hypix_θΨK,
					Plot_HypixTime,
					Plot_RainfallInterception,
					Plot_Se_Time,
					Plot_Se_Z,
					Plot_Sorptivity,
					Vegetation,

					Plot_OfStep,
					Plot_Se_Ψ_Constrained,
					Plot_θΨ_Δθ,
					Plot_σ2θr,
					Plot_Ψmin_Ψmax,
					Plot_θ∂θ∂Ψ)

		path = PATHS(Home, hyPix, inputSmap, inputSoilwater, inputTemporary, option, plotSoilwater, smap2Hypix, tableSmap, tableSoilwater)

	return path
	end # function PATHS			
end  # module path
# ............................................................