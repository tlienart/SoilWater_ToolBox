# =============================================================
#		MODULE: path
# =============================================================
module path

	# INPUT PATH
		Hkunsat = "H_Kunsat.csv"
		Hθ = "H_Theta.csv"
		Psd = "Psd.csv"
		Infiltration = "Infiltration.csv"
		IdTrue = "Id_True.csv"
		Parameter = "Parameter.csv"

	# OUTPUT PATH

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# PROCESSING 
		Home = @__DIR__

		Hkunsat = Home * "//INPUT//" * Hkunsat
		Hθ = Home * "//INPUT//" * Hθ
		Psd = Home * "//INPUT//" * Psd
		Infiltration = Home * "//INPUT//" * Infiltration
		IdTrue = Home * "//INPUT//" * IdTrue
end  # module path
# ............................................................
