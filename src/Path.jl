# =============================================================
#		MODULE: path
# =============================================================
module paths

	mutable struct INPUT_SMAP
		HydroParam_ThetaH::String
		Smap::String
		SmapLookupTableWettable::String
	end

	mutable struct INPUT_SOILWATER
		BulkDensity::String
		ConvertModel::String
		HydroParam_Infilt::String
		HydroParam_ThetaH::String 
		Id_Select::String
		Infiltration_Param::String
		Infiltration::String
		Kunsat_Model::String
		Kunsat::String
		Psd::String
		Ψθ::String
	end # struct INPUT_SOILWATER

	mutable struct INPUT_TEMPORARY
		σ_ψM_Scenario::String
	end
	
	struct OPTION 
		Model_Name::String
		ProjectName_Hypix::String
		Select::String
		SiteName_Soilhyro::String
	end # struct OPTION
	
	mutable struct PLOT_SOILWATER	
		Plot_Psd_θΨ::String				
		Plot_∑infilt_Opt::String
		Plot_∑infilt_SeIniRange::String
		Plot_∑infilt_θΨ::String
		Plot_IMP_model::String
		Plot_Psd_θr::String
		Plot_Sorptivity_Se::String
		Plot_θΨK::String
		Plot_σΨm::String
	end # PLOT_SOILWATER

	
	mutable struct TABLE_SMAP
		Table_θΨK₀           = "Table_ThetaHK.csv"
	end # struct INPUT_TABLE_SMAP

	mutable struct TABLE_SOILWATER
		Table_ExtraPoints_θΨ::String
		Table_HydroInfilt::String
		Table_Infilt::String
		Table_KosugiθΨ::String
		Table_Psd_θΨ_θ::String
		Table_Psd::String
		Table_θΨ_Psd::String
	end # struct INPUT_SOILWATER

	mutable struct PATHS
		inputSmap::INPUT_SMAP
		inputSoilwater::INPUT_SOILWATER
		inputTemporary::INPUT_TEMPORARY
		option::OPTION
		plotSoilwater::PLOT_SOILWATER
		tableSmap::TABLE_SMAP
		tableSoilwater::INPUT_SOILWATER
	end


	#__________________________________________________________________
	#..................................................................

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PATHS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PATHS()

		# Change path name to /data/Private/
			Home2 = @__DIR__

		# perform cs..
			Home = dirname(Home2)

			Home = Home * "/data/Private/"

		# Paths
			FileDataSoilhydro_Input = Home * "INPUT//DataSoilHydraulic//" * SiteName_Soilhyro * "//" * SiteName_Soilhyro * "_"

			FileSoilHydro_Table₁ = Home * "//OUTPUT//SoilHydro//" * SiteName_Soilhyro * "//Table//" 
				mkpath(FileSoilHydro_Table₁) 


		# =============================================================
		#		OPTIONS
		# =============================================================
			# Which files to use
				SiteName_Soilhyro = "Smap20210226" #"VCSNSmap2"; "SFF"; "PAF"; K10KPA; Smap; Smap20210226; SmapSouthland2; CantyLysimSmap; VCSNSmap; "WaikLysim"; "Convert; "SmapNZAllSoilsSmap20210326"; "Smap20210226"
				Model_Name ="A"
				Select = "SELECT_1" # Select data to model

			option = OPTION(Model_Name, ProjectName_Hypix, Select, SiteName_Soilhyro)

		# =============================================================
		#		INPUT_SOILWATER
		# =============================================================
			# DATA input into SoilWater
				BulkDensity             = "BulkDensity.csv"
				ConvertModel            = "TableHydro_Compiled_Homogeneous.csv"
				HydroParam_Infilt       = "HydroParam_Infilt.csv"
				HydroParam_ThetaH       = "GUI_HydroParam.csv"
				Id_Select               = "IdSelect.csv"
				Infiltration            = "Infiltration.csv"
				Infiltration_Param      = "Infiltration_Param.csv"
				Kunsat                  = "KunsatH.csv"
				Kunsat_Model            = "Kunsat_H_model.csv"
				Psd                     = "Psd.csv"
				Ψθ                      = "ThetaH.csv"

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            BulkDensity        = FileDataSoilhydro_Input * BulkDensity
            ConvertModel       = FileDataSoilhydro_Input * ConvertModel
            HydroParam_Infilt  = FileDataSoilhydro_Input * HydroParam_Infilt
            HydroParam_ThetaH  = FileDataSoilhydro_Input * HydroParam_ThetaH
            Id_Select          = FileDataSoilhydro_Input * Id_Select
            Infiltration       = FileDataSoilhydro_Input * Infiltration
            Infiltration_Param = FileDataSoilhydro_Input * Infiltration_Param
            Kunsat             = FileDataSoilhydro_Input * Kunsat
            Kunsat_Model       = FileDataSoilhydro_Input * Kunsat_Model
            Psd                = FileDataSoilhydro_Input * Psd
            Ψθ                 = FileDataSoilhydro_Input * Ψθ

			inputSoilwater = INPUT_SOILWATER(BulkDensity, ConvertModel, HydroParam_Infilt, HydroParam_ThetaH, Id_Select, Infiltration, Infiltration_Param, Kunsat, Kunsat_Model, Psd, Ψθ)		

				
		# =============================================================
		#		INPUT_SMAP
		#		None core
		# =============================================================
			# Smap input path				
				Smap                    = "Layer.csv"
				SmapLookupTableWettable = "LookupTable_Stone.csv"

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				Smap = FileDataSoilhydro_Input * Smap
				SmapLookupTableWettable = FileDataSoilhydro_Input * SmapLookupTableWettable

			inputSmap = INPUT_SMAP(HydroParam_ThetaH, Smap, SmapLookupTableWettable)
				

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
				Table_ExtraPoints_θΨ = "Table_ExtraPoints_θΨ.csv"
				Table_HydroInfilt    = "Table_HydroInfilt.csv"
				Table_Infilt         = "Table_Infilt.csv"
				Table_KosugiθΨ       = "Table_KosugiθΨ.csv"
				Table_Psd            = "Table_Psd.csv"
				Table_Psd_θΨ_θ       = "Table_PsdTheta.csv"
				Table_θΨ_Psd         = "Table_PsdHydro.csv"
			
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				FileSoilHydro_Table₁ = FileSoilHydro_Table₁ * SiteName_Soilhyro * "_"

				Table_ExtraPoints_θΨ = FileSoilHydro_Table₁ *   "_" * Table_ExtraPoints_θΨ
				Table_HydroInfilt    = FileSoilHydro_Table₁ * string(option.infilt.Model) * "_" *  Model_Name  *  "_" * Table_HydroInfilt
				Table_Infilt         = FileSoilHydro_Table₁ * string(option.infilt.Model) *  "_" *  Model_Name  *  "_" *  Table_Infilt
				Table_KosugiθΨ       = FileSoilHydro_Table₁ *   "_" * Table_KosugiθΨ
				Table_Psd            = FileSoilHydro_Table₁ * string(option.psd.Model) *  "_" * Model_Name * "_" * Table_Psd
				Table_Psd_θΨ_θ       = FileSoilHydro_Table₁ * string(option.psd.HydroModel) *  "_" * Model_Name * "_" *  Table_Psd_θΨ_θ
				Table_θΨ_Psd         = FileSoilHydro_Table₁ * string(option.psd.HydroModel) *  "_" * string(option.hydro.σ_2_Ψm) *  "_" * Model_Name * "_" * Table_θΨ_Psd
			
			tableSoilwater = TABLE_SOILWATER(Table_ExtraPoints_θΨ, Table_HydroInfilt, Table_Infilt, Table_KosugiθΨ, Table_Psd, Table_Psd_θΨ_θ, Table_θΨ_Psd)


		# =============================================================
		#		TABLE_SMAP
		# =============================================================
			# Output tables smap
				Table_θΨK₀           = "Table_ThetaHK.csv"

				Table_θΨK            = FileSoilHydro_Table₁ * string(option.hydro.HydroModel) *  "_" * string(option.hydro.σ_2_Ψm) *  "_" * Model_Name * "_" * Table_θΨK₀

				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				Table_Smap =  FileSoilHydro_Table₁ * "Smap.csv"

				tableSmap = TABLE_SMAP(Table_θΨK)

				
		# =============================================================
		#		PLOT SOILWATER
		# =============================================================
			FileSoilHydro_Plot = Home * "//OUTPUT//SoilHydro//" * SiteName_Soilhyro * "//Plots//"

			Plot_θΨK = FileSoilHydro_Plot * "//Lab//" 
				mkpath(Plot_θΨK)
				Plot_θΨK = Plot_θΨK * SiteName_Soilhyro * "_"

			Plot_σΨm = FileSoilHydro_Plot * "//LabSigmaHm//" 
				mkpath(Plot_σΨm)
				Plot_σΨm = Plot_σΨm * SiteName_Soilhyro * "_"

			Plot_Psd_θΨ = FileSoilHydro_Plot * "//Psd//IMP_ThetaH//"
				mkpath(Plot_Psd_θΨ)				
				Plot_Psd_θΨ = Plot_Psd_θΨ * SiteName_Soilhyro * "_"

			Plot_IMP_model = FileSoilHydro_Plot * "//Psd//IMP//"
				mkpath(Plot_IMP_model)
				Plot_IMP_model = Plot_IMP_model * SiteName_Soilhyro * "_"

			Plot_Psd_θr = FileSoilHydro_Plot * "//Psd//ThetaR//" 
				mkpath(Plot_Psd_θr)
				Plot_Psd_θr = Plot_Psd_θr * "Plot_ThetaR.svg"

			Plot_∑infilt_Opt = FileSoilHydro_Plot * "//Infiltration//Optimize//"
				mkpath(Plot_∑infilt_Opt)
				Plot_∑infilt_Opt = Plot_∑infilt_Opt * SiteName_Soilhyro * "_"

			Plot_∑infilt_SeIniRange = FileSoilHydro_Plot * "//Infiltration//SeIni//"
				mkpath(Plot_∑infilt_SeIniRange)
				Plot_∑infilt_SeIniRange = Plot_∑infilt_SeIniRange * SiteName_Soilhyro * "_"

			Plot_∑infilt_θΨ = FileSoilHydro_Plot * "//Infiltration//ThetaH//"
				mkpath(Plot_∑infilt_θΨ)
				Plot_∑infilt_θΨ = Plot_∑infilt_θΨ * SiteName_Soilhyro * "_"

			Plot_Sorptivity_Se = FileSoilHydro_Plot * "//Infiltration//Sorptivity//"
				mkpath(Plot_Sorptivity_Se)
				Plot_Sorptivity_Se = Plot_Sorptivity_Se * SiteName_Soilhyro * "_"

		plotSoilwater = PLOT_SOILWATER(Plot_Psd_θΨ, Plot_∑infilt_Opt, Plot_∑infilt_SeIniRange, Plot_∑infilt_θΨ, Plot_IMP_model, Plot_Psd_θr, Plot_Sorptivity_Se, Plot_θΨK, Plot_σΨm)	

		path = PATHS(inputSmap, inputSoilwater, inputTemporary, option, plotSoilwater, tableSmap, tableSoilwater)

	return path
	end # function PATHS			
end  # module path
# ............................................................