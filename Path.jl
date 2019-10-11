# =============================================================
#		MODULE: path
# =============================================================
module path

	# NAME OF FILE
		Name = "PAF_"

	# INPUT PATH
		Id_Select 			= "Id_Select.csv"
		Kunsat 				= "Kunsat_H.csv"
		Ψθ 					= "Theta_H.csv"
		Psd 				= "Psd.csv"
		Infiltration 		= "Infiltration_1.csv"
		Infiltration_Param 	= "Infiltration_Param.csv"


	# OUTPUT PATH
		Table_θΨK			= "Table_ThetaHK.csv"
		
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# PROCESSING 
		Home = @__DIR__
		Id_Select = Home * "//INPUT//" * Name * Id_Select
		Kunsat = Home * "//INPUT//" * Name * Kunsat
		Ψθ = Home * "//INPUT//" * Name * Ψθ
		Psd = Home * "//INPUT//"  * Name * Psd
		Infiltration = Home * "//INPUT//" * Name * Infiltration
		Infiltration_Param = Home * "//INPUT//" * Name * Infiltration_Param

		Table_θΨK = Home * "//OUTPUT//Table//" * Name * Table_θΨK
		Plots_θΨK = Home * "//OUTPUT//Plots//"  * Name
end  # module path
# ............................................................
