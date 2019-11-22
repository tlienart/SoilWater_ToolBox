# =============================================================
#		MODULE: path
# =============================================================
module path
	import ..option

	# NAME OF FILE
		Name = "SFF_"

	# INPUT PATH
		Id_Select 			= "Id_Select.csv"
		Kunsat 				= "Kunsat_H.csv"
		Ψθ 					= "Theta_H.csv"
		Psd 				= "Psd.csv"
		Infiltration 		= "Infiltration_1.csv"
		Infiltration_Param 	= "Infiltration_Param.csv"
		PsdΦ		 		= "PsdPorosity.csv" 

	# OUTPUT PATH
		Table_θΨK			= "Table_ThetaHK.csv"
		Table_Psd			= "Table_Psd.csv"
		Table_θΨK_Psd		= "Table_PsdHydro.csv"
		Table_θr			= "Table_Thetar.csv"
		
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# PROCESSING 
	
		Home = @__DIR__
		Id_Select = Home * "//INPUT//" * Name * Id_Select
		Kunsat = Home * "//INPUT//" * Name * Kunsat
		Ψθ = Home * "//INPUT//" * Name * Ψθ
		Psd = Home * "//INPUT//"  * Name * Psd
		Infiltration = Home * "//INPUT//" * Name * Infiltration
		Infiltration_Param = Home * "//INPUT//" * Name * Infiltration_Param
		PsdΦ = Home * "//INPUT//" * Name * PsdΦ

        Table_θΨK                = Home * "//OUTPUT//Table//" * Name *  option.hydro.HydroModel * "_" * Table_θΨK
        Table_Psd                = Home * "//OUTPUT//Table//" * Name * option.psd.Model * "_" * Table_Psd
        Table_θr                 = Home * "//OUTPUT//Table//" * Name * Table_θr
        Table_θΨK_Psd            = Home * "//OUTPUT//Table//" * Name * option.psd.HydroModel * "_" * Table_θΨK_Psd
        Plots_θΨK                = Home * "//OUTPUT//Plots//Lab//" * Name
        Plots_BestLab            = Home * "//OUTPUT//Plots//Infiltration//Lab//" * Name
        Plots_BestLab_SeIniRange = Home * "//OUTPUT//Plots//Infiltration//Lab_SeIniRange//" * Name
        Plots_Psd                = Home * "//OUTPUT//Plots//Psd//" * Name
        Plots_Psd_ThetaR         = Home * "//OUTPUT//Plots//Psd//ThetaR//" * Name * "_Plot_ThetaR.svg"
        Plots_IMP_model          = Home * "//OUTPUT//Plots//Psd//IMP_results//" * Name


end  # module path
# ............................................................
