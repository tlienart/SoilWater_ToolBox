# =============================================================
#		MODULE: path
# =============================================================
module path

	# INPUT PATH
		Id_Select 			= "Id_Select.csv"
		Kunsat 				= "Kunsat_H.csv"
		Ψθ 					= "Theta_H.csv"
		Psd 				= "Psd.csv"
		Infiltration 		= "Infiltration.csv"
		Infiltration_Param 	= "Infiltration_Param.csv"
		Parameter 			= "Parameter.csv"

	# OUTPUT PATH
		Table_θΨK			= "Table_ThetaHK.csv"

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# PROCESSING 
		Home = @__DIR__
		Id_Select = Home * "//INPUT//" * Id_Select
		Kunsat = Home * "//INPUT//" * Kunsat
		Ψθ = Home * "//INPUT//" * Ψθ
		Psd = Home * "//INPUT//" * Psd
		Infiltration = Home * "//INPUT//" * Infiltration
		Infiltration_Param = Home * "//INPUT//" * Infiltration_Param

		Table_θΨK = Home * "//OUTPUT//Table//" * Table_θΨK
end  # module path
# ............................................................
