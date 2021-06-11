# =============================================================
#		MODULE: path
# =============================================================
module paths

	struct OPTION 
		SiteName_Soilhyro::String
		ProjectName_Hypix::String
		Model_Name::String
		Select::String
	end # struct OPTION 

	mutable struct INPUT_SOILWATER
         ConvertModel::String
         Id_Select::String
         Infiltration::String
         Infiltration_Param::String
         Kunsat::String
         Kunsat_Model::String
         Psd::String
         Ψθ::String
         BulkDensity::String
         HydroParam_Infilt::String
	end # struct INPUT_SOILWATER

	mutable struct INPUT_SMAP
		Smap::String
		HydroParam_ThetaH::String
		SmapLookupTableWettable::String
	end

	mutable struct INPUT_TEMPORARY
		σ_ψM_Scenario::String
	end
	
	mutable struct INPUT_HYPIX
		 a
	end # struct INPUT_SOILWATER
	
	mutable struct TABLE_SOILWATER
		Table_HydroInfilt::String
		Table_Infilt::String
		Table_Psd::String
		Table_Psd_θΨ_θ::String
		Table_θΨ_Psd::String
		Table_ExtraPoints_θΨ::String
		Table_KosugiθΨ::String
	end # struct INPUT_SOILWATER

	mutable struct PLOT_SOILWATER	
		Plots_θΨK::String
		Plots_σΨm::String
		Plot_Psd_θΨ::String				
		Plots_IMP_model::String
		Plots_Psd::String
		Plots_Psd_θr::String
		Plots_∑infilt_Opt::String
		Plots_∑infilt_SeIniRange::String
		Plots_∑infilt_θΨ::String
		Plots_Sorptivity_Se::String
	end # PLOT_SOILWATER

	mutable struct TABLE_SMAP
		Table_θΨK₀           = "Table_ThetaHK.csv"
	end # struct INPUT_TABLE_SMAP
	
	mutable struct OUTPUT_HYPIX
		 
	end # struct INPUT_SOILWATER

	mutable struct PATHS
		option::OPTION
		inputSoilwater::INPUT_SOILWATER
		inputSmap::INPUT_SMAP
		inputTemperory::INPUT_TEMPORARY
		inputHypix::INPUT_HYPIX

		tableSoilwater::INPUT_SOILWATER
		tableSmap::INPUT_SMAP
		tableTemperory::INPUT_TEMPORARY
		tableHypix::INPUT_HYPIX
		plotSoilwater::PLOT_SOILWATER
	end

	#__________________________________________________________________
	#..................................................................

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PATHS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PATHS()

		# NAME OF FILE
			SiteName_Soilhyro = "Smap20210226" #"VCSNSmap2"; "SFF"; "PAF"; K10KPA; Smap; Smap20210226; SmapSouthland2; CantyLysimSmap; VCSNSmap; "WaikLysim"; "Convert; "SmapNZAllSoilsSmap20210326"; "Smap20210226"
			ProjectName_Hypix = "JULES" # "JULES"; "LYSIMETERS" 
			Model_Name ="A"
			Select = "SELECT_1" # Select data to model
			
		# INPUT PATH
			# Smap
				Smap                    = "Layer.csv"
				HydroParam_ThetaH       = "GUI_HydroParam.csv"

			# DATA SoilWater_Toolbox
			ConvertModel            = "TableHydro_Compiled_Homogeneous.csv"
			SmapLookupTableWettable = "LookupTable_Stone.csv"
			Id_Select               = "IdSelect.csv"
			Infiltration            = "Infiltration.csv"
			Infiltration_Param      = "Infiltration_Param.csv"
			Kunsat                  = "KunsatH.csv"
			Kunsat_Model            = "Kunsat_H_model.csv"
			Psd                     = "Psd.csv"
			PsdΦ                    = "PsdPorosity.csv"
			Ψθ                      = "ThetaH.csv"
			BulkDensity             = "BulkDensity.csv"
			HydroParam_Infilt       = "HydroParam_Infilt.csv"
			σ_ψM_Scenario           = "σ_ψM_Scenario.csv"


		# TABLE OUTPUT PATH
			# DATA SOIL HYDRO
			Table_HydroInfilt    = "Table_HydroInfilt.csv"
			Table_Infilt         = "Table_Infilt.csv"
			Table_Psd            = "Table_Psd.csv"
			Table_Psd_θΨ_θ       = "Table_PsdTheta.csv"
			Table_θΨK₀           = "Table_ThetaHK.csv"
			Table_θΨ_Psd         = "Table_PsdHydro.csv"
			Table_ExtraPoints_θΨ = "Table_ExtraPoints_θΨ.csv"
			Table_KosugiθΨ       = "Table_KosugiθΨ.csv"

		#__________________________________________________________________
		#..................................................................
		#__________________________________________________________________

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PATHS
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PATHS()
		
		# PROCESSING
			# Change path name to /data/Private/
			Home2 = @__DIR__

			# perform cs..
			Home = dirname(Home2)

			Home = Home * "/data/Private/"

			# INPUT
				# DATA SOIL HYDRO
				FileDataSoilhydro_Input = Home * "INPUT//DataSoilHydraulic//" * SiteName_Soilhyro * "//" * SiteName_Soilhyro * "_"

				# ID_Select	
					Id_Select          = FileDataSoilhydro_Input * Id_Select

					# Smap
						Smap = FileDataSoilhydro_Input * Smap
						SmapLookupTableWettable = FileDataSoilhydro_Input * SmapLookupTableWettable

					# Convert
						ConvertModel = FileDataSoilhydro_Input * ConvertModel
						
					# BulkDensity
						BulkDensity       = FileDataSoilhydro_Input * BulkDensity
					
					#Lab
				Ψθ                = FileDataSoilhydro_Input * Ψθ
				Kunsat            = FileDataSoilhydro_Input * Kunsat
				Kunsat_Model      = FileDataSoilhydro_Input * Kunsat_Model
				HydroParam_ThetaH = FileDataSoilhydro_Input * HydroParam_ThetaH
				σ_ψM_Scenario     = FileDataSoilhydro_Input * σ_ψM_Scenario
				Temporary_1       = FileDataSoilhydro_Input * Temporary_1
				Temporary_2       = FileDataSoilhydro_Input * Temporary_2
						
					#Psd
				Psd            = FileDataSoilhydro_Input * Psd
				# HydroParam_Psd = FileDataSoilhydro_Input * HydroParam_Psd

					# Infilt
				Infiltration       = FileDataSoilhydro_Input * Infiltration
				Infiltration_Param = FileDataSoilhydro_Input * Infiltration_Param
				HydroParam_Infilt  = FileDataSoilhydro_Input * HydroParam_Infilt

		
			# TABLE
				# SOIL HYDRO
				FileSoilHydro_Table₁ = Home * "//OUTPUT//SoilHydro//" * SiteName_Soilhyro * "//Table//" 
				#Make Folder if not exist
				mkpath(FileSoilHydro_Table₁) 
				FileSoilHydro_Table₁ = FileSoilHydro_Table₁ * SiteName_Soilhyro * "_"

					#Lab
				Table_θΨK            = FileSoilHydro_Table₁ * string(option.hydro.HydroModel) *  "_" * string(option.hydro.σ_2_Ψm) *  "_" * Model_Name * "_" * Table_θΨK₀
				Table_θΨ_Psd         = FileSoilHydro_Table₁ * string(option.psd.HydroModel) *  "_" * string(option.hydro.σ_2_Ψm) *  "_" * Model_Name * "_" * Table_θΨ_Psd
				Table_ExtraPoints_θΨ = FileSoilHydro_Table₁ *   "_" * Table_ExtraPoints_θΨ
				Table_KosugiθΨ       = FileSoilHydro_Table₁ *   "_" * Table_KosugiθΨ

					#SMAP
						Table_Smap =  FileSoilHydro_Table₁ * "Smap.csv"
	
					#Infilt
				Table_HydroInfilt    = FileSoilHydro_Table₁ * string(option.infilt.Model) * "_" *  Model_Name  *  "_" * Table_HydroInfilt
				Table_Infilt         = FileSoilHydro_Table₁ * string(option.infilt.Model) *  "_" *  Model_Name  *  "_" *  Table_Infilt

					#Psd
						Table_Psd         = FileSoilHydro_Table₁ * string(option.psd.Model) *  "_" * Model_Name * "_" * Table_Psd
						Table_Psd_θΨ_θ    = FileSoilHydro_Table₁ * string(option.psd.HydroModel) *  "_" * Model_Name * "_" *  Table_Psd_θΨ_θ

				
				
			# PLOT
				# SOIL HYDRO
				FileSoilHydro_Plot = Home * "//OUTPUT//SoilHydro//" * SiteName_Soilhyro * "//Plots//"
					#Lab
				Plots_θΨK = FileSoilHydro_Plot * "//Lab//" 
						mkpath(Plots_θΨK)
						Plots_θΨK  = Plots_θΨK * SiteName_Soilhyro * "_"

						Plots_σΨm = FileSoilHydro_Plot * "//LabSigmaHm//" 
						mkpath(Plots_σΨm)
						Plots_σΨm  = Plots_σΨm * SiteName_Soilhyro * "_"


					#Psd
				Plot_Psd_θΨ     = FileSoilHydro_Plot * "//Psd//IMP_ThetaH//"
						mkpath(Plot_Psd_θΨ)				
				Plot_Psd_θΨ = Plot_Psd_θΨ * SiteName_Soilhyro * "_"
						
				Plots_IMP_model = FileSoilHydro_Plot * "//Psd//IMP//"
						mkpath(Plots_IMP_model)
				Plots_IMP_model = Plots_IMP_model * SiteName_Soilhyro * "_"
						
				Plots_Psd       = FileSoilHydro_Plot * "//Psd//"
						mkpath(Plots_Psd)
				Plots_Psd       = Plots_Psd * SiteName_Soilhyro * "_"

				Plots_Psd_θr    = FileSoilHydro_Plot * "//Psd//ThetaR//" 
						mkpath(Plots_Psd_θr)
						Plots_Psd_θr    = Plots_Psd_θr * "Plot_ThetaR.svg"
						
					#Infiltration					
				Plots_∑infilt_Opt        = FileSoilHydro_Plot * "//Infiltration//Optimize//"
						mkpath(Plots_∑infilt_Opt)
				Plots_∑infilt_Opt        = Plots_∑infilt_Opt * SiteName_Soilhyro * "_"
						
				Plots_∑infilt_SeIniRange = FileSoilHydro_Plot * "//Infiltration//SeIni//"
						mkpath(Plots_∑infilt_SeIniRange)
				Plots_∑infilt_SeIniRange = Plots_∑infilt_SeIniRange * SiteName_Soilhyro * "_"

				Plots_∑infilt_θΨ         = FileSoilHydro_Plot * "//Infiltration//ThetaH//"
						mkpath(Plots_∑infilt_θΨ)
				Plots_∑infilt_θΨ         = Plots_∑infilt_θΨ * SiteName_Soilhyro * "_"
						
				Plots_Sorptivity_Se      = FileSoilHydro_Plot * "//Infiltration//Sorptivity//"
						mkpath(Plots_Sorptivity_Se)
				Plots_Sorptivity_Se      = Plots_Sorptivity_Se * SiteName_Soilhyro * "_"
			
			paths = PATHS(option, inputSoilwater, inputHypix, outputSoilwater, outputHypix)
	return path
	end # function PATHS			
end  # module path
# ............................................................
