# =============================================================
#		MODULE: path
# =============================================================
module path

	# INPUT PATH
		Id_Select 		= "Id_Select.csv"
		Kunsat 			= "Kunsat_H.csv"
		Ψθ 				= "Theta_H.csv"
		Psd 			= "Psd.csv"
		Infiltration 	= "Infiltration.csv"
		Parameter 		= "Parameter.csv"

	# OUTPUT PATH

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# PROCESSING 
		Home = @__DIR__
		Id_Select = Home * "//INPUT//" * Id_Select

		Kunsat = Home * "//INPUT//" * Kunsat
		Ψθ = Home * "//INPUT//" * Ψθ
		Psd = Home * "//INPUT//" * Psd
		Infiltration = Home * "//INPUT//" * Infiltration
end  # module path
# ............................................................
