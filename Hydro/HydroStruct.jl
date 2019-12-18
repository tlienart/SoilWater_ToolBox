# =============================================================
#		MODULE: hydroStruct
# =============================================================
module hydroStruct
	import ..option, ..tool

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : HYDRAULIC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mutable struct KOSUGI
            θs     :: 	Vector{Float64}
            θr     :: 	Vector{Float64}
            σ      :: 	Vector{Float64}
            Ψm     :: 	Vector{Float64}
            Ks     :: 	Vector{Float64}
            θsMat  ::	Vector{Float64}
            σMac   :: 	Vector{Float64}
				ΨmMac  ::	Vector{Float64}
				Φ  	 ::	Vector{Float64}
            Nse    :: 	Vector{Float64}
            Nse_θψ :: 	Vector{Float64}
            Nse_Kψ :: 	Vector{Float64}

				FieldName ::	Vector{Symbol} # Need to put
		end # struct KOSUGI

		mutable struct VANGENUCHTEN
         θs        ::	Vector{Float64}
         θr        ::	Vector{Float64}
         N         ::	Vector{Float64}
         Ψvg       ::	Vector{Float64}
         Ks        ::	Vector{Float64}
         Φ         ::	Vector{Float64}
         Nse       ::	Vector{Float64}
         Nse_θψ    ::	Vector{Float64}
         Nse_Kψ    ::	Vector{Float64}

         FieldName ::	Vector{Symbol} # Need to put
		end # struct VANGENUCHTEN


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROSTRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDROSTRUCT(N_SoilSelect)
			# For all models
			FieldName = Array{Symbol}(undef, 1) # Need to put
			θs        = zeros(Float64, N_SoilSelect)
			θr        = zeros(Float64, N_SoilSelect)
			Ks        = zeros(Float64, N_SoilSelect)
			θsMat     = zeros(Float64, N_SoilSelect)
			Φ    	  = zeros(Float64, N_SoilSelect)
			Nse       = zeros(Float64, N_SoilSelect)
			Nse_θψ    = zeros(Float64, N_SoilSelect)
			Nse_Kψ    = zeros(Float64, N_SoilSelect)

			if option.hydro.HydroModel == "Kosugi" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
                σ     = zeros(Float64, N_SoilSelect)
                Ψm    = zeros(Float64, N_SoilSelect)
                σMac  = zeros(Float64, N_SoilSelect)
                ΨmMac = zeros(Float64, N_SoilSelect)

				hydro = KOSUGI(θs, θr, σ, Ψm, Ks, θsMat, σMac, ΨmMac, Φ, Nse, Nse_θψ, Nse_Kψ, FieldName)

				return hydro = tool.readWrite.FIELDNAME_2_STRUCT(KOSUGI, hydro) # Saving the FieldNames

	
			elseif option.hydro.HydroModel == "Vangenuchten" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<><>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				N		= zeros(Float64, N_SoilSelect)
				Ψvg		= zeros(Float64, N_SoilSelect)

				hydro = VANGENUCHTEN(θs, θr, N, Ψvg, Ks, Φ, Nse, Nse_θψ, Nse_Kψ, FieldName) # Need to put

				return hydro = tool.readWrite.FIELDNAME_2_STRUCT(VANGENUCHTEN, hydro) # Saving the FieldNames
			end # option.hydro.HydroModel

		end #  function HYDROSTRUCT
end  # module: hydroStruct

# ............................................................