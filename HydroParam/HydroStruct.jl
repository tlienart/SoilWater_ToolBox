# =============================================================
#		MODULE: hydroStruct
# =============================================================
module hydroStruct
	import ..option

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : HYDRAULIC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mutable struct KOSUGI
			θs :: 		Vector{Float64}
			θr :: 		Vector{Float64}
			Ks :: 		Vector{Float64}
			σ :: 		Vector{Float64}
			Ψm :: 		Vector{Float64}
			θsMat ::	Vector{Float64}
			σMac :: 	Vector{Float64}
			ΨmMac ::	Vector{Float64}
		end # struct KOSUGI

		mutable struct VANGENUCHTEN
			θs :: 	Vector{Float64}
			θr :: 	Vector{Float64}
			Ks :: 	Vector{Float64}
			N :: 	Vector{Float64}
			Ψvg :: 	Vector{Float64}
		end # struct VANGENUCHTEN


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROSTRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
		function HYDROSTRUCT(N_SoilSelect)
			# For all models
				θs 		= Array{Float64}(undef, (N_SoilSelect))
				θr 		= Array{Float64}(undef, (N_SoilSelect))
				Ks		= zeros(Float64, N_SoilSelect)
				θsMat 	= Array{Float64}(undef, (N_SoilSelect))

			if option.HydroModel == "Kosugi"
				σ 		= Array{Float64}(undef, (N_SoilSelect))
				Ψm 		= Array{Float64}(undef, (N_SoilSelect))
				σMac 	= Array{Float64}(undef, (N_SoilSelect))
				ΨmMac 	= Array{Float64}(undef, (N_SoilSelect))
				
				return hydro = KOSUGI(θs, θr, Ks, σ, Ψm, θsMat, σMac, ΨmMac)
				
			elseif option.HydroModel == "Vangenuchten"
				N		= Array{Float64}(undef, (N_SoilSelect))
				Ψvg		= Array{Float64}(undef, (N_SoilSelect))

				return hydro = VANGENUCHTEN(θs, θr, Ks, N, Ψvg)
			end # option.HydroModel

		end #  function HYDROSTRUCT
end  # module: hydroStruct

# ............................................................
