# =============================================================
#		MODULE: path
# =============================================================
module path
	import ..option

	# NAME OF FILE
		Name = "PAF_"

	# INPUT PATH
      Id_Select          = "Id_Select.csv"
      Infiltration       = "Infiltration_1.csv"
      Infiltration_Param = "Infiltration_Param.csv"
      Kunsat             = "Kunsat_H.csv"
      Psd                = "Psd.csv"
      PsdΦ               = "PsdPorosity.csv"
      Ψθ                 = "Theta_H.csv"
      ρ_Ψθ               = "BulkDensity_ThetaH.csv"
      ρ_Psd              = "BulkDensity_Psd.csv"
      ρ_Infilt           = "BulkDensity_Infilt.csv"

	# OUTPUT PATH
      Table_HydroInfilt = "Table_HydroInfilt.csv"
      Table_Infilt      = "Table_Infilt.csv"
      Table_Psd         = "Table_Psd.csv"
      Table_θr          = "Table_Thetar.csv"
      Table_θΨK         = "Table_ThetaHK.csv"
      Table_θΨK_Psd     = "Table_PsdHydro.csv"
		
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# PROCESSING
		Home = @__DIR__

		# Input
         Id_Select          = Home * "//INPUT//" * Name * Id_Select
         Infiltration       = Home * "//INPUT//" * Name * Infiltration
         Infiltration_Param = Home * "//INPUT//" * Name * Infiltration_Param
         Kunsat             = Home * "//INPUT//" * Name * Kunsat
         Psd                = Home * "//INPUT//"  * Name * Psd
         PsdΦ               = Home * "//INPUT//" * Name * PsdΦ
         Ψθ                 = Home * "//INPUT//" * Name * Ψθ
         ρ_Infilt           = Home * "//INPUT//" * Name * ρ_Infilt
         ρ_Psd              = Home * "//INPUT//" * Name * ρ_Psd
         ρ_Ψθ               = Home * "//INPUT//" * Name * ρ_Ψθ

		# Table
			Table_θΨK                = Home * "//OUTPUT//Table//" * Name *  option.hydro.HydroModel * "_" * Table_θΨK
			Table_Psd                = Home * "//OUTPUT//Table//" * Name * option.psd.Model * "_" * Table_Psd
			Table_θr                 = Home * "//OUTPUT//Table//" * Name * Table_θr
			Table_θΨK_Psd            = Home * "//OUTPUT//Table//" * Name * option.psd.HydroModel * "_" * Table_θΨK_Psd
			Table_HydroInfilt        = Home * "//OUTPUT//Table//" * Name * Table_HydroInfilt
			Table_Infilt             = Home * "//OUTPUT//Table//" * Name * option.infilt.Model * "_" *  Table_Infilt
			Plots_θΨK                = Home * "//OUTPUT//Plots//Lab//" * Name
			Plots_BestLab            = Home * "//OUTPUT//Plots//Infiltration//Lab//" * Name
			Plots_BestLab_SeIniRange = Home * "//OUTPUT//Plots//Infiltration//Lab_SeIniRange//" * Name
			Plots_∑infilt_Tinfilt    = Home * "//OUTPUT//Plots//Infiltration//" * Name
			Plots_Psd                = Home * "//OUTPUT//Plots//Psd//" * Name
			Plots_Psd_θr         = Home * "//OUTPUT//Plots//Psd//ThetaR//" * Name * "_Plot_ThetaR.svg"
			Plots_IMP_model          = Home * "//OUTPUT//Plots//Psd//IMP_results//" * Name


end  # module path
# ............................................................
