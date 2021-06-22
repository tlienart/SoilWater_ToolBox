# =============================================================
#		MODULE: hydroStruct
# =============================================================
module hydroStruct
	import ..tool
	export HYDROSTRUCT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : HYDRAULIC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mutable struct KOSUGI  # <>=<>=<>=<>=<>=<>=<>=<>=<>
         θs             :: Vector{Float64}
         θr             :: Vector{Float64}
         Ks             :: Vector{Float64}
         Ψm             :: Vector{Float64}
         σ              :: Vector{Float64}
         θsMacMat_ƞ     ::	Vector{Float64}
         σMac           :: Vector{Float64}
         ΨmMac          ::	Vector{Float64}
         So             :: Vector{Float64}
         θsMacMat       ::	Vector{Float64}
         Φ              ::	Vector{Float64}
			
         θs_Min         :: Vector{Float64}
         θr_Min         :: Vector{Float64}
         Ks_Min         :: Vector{Float64}
         Ψm_Min         :: Vector{Float64}
         σ_Min          :: Vector{Float64}
         θsMacMat_ƞ_Min ::	Vector{Float64}
         σMac_Min       :: Vector{Float64}
         ΨmMac_Min      ::	Vector{Float64}
         So_Min         :: Vector{Float64}
         θsMacMat_Min   ::	Vector{Float64}
         Φ_Min          ::	Vector{Float64}

         θs_Max         :: Vector{Float64}
         θr_Max         :: Vector{Float64}
         Ks_Max         :: Vector{Float64}
         Ψm_Max         :: Vector{Float64}
         σ_Max          :: Vector{Float64}
         θsMacMat_ƞ_Max ::	Vector{Float64}
         σMac_Max       :: Vector{Float64}
         ΨmMac_Max      ::	Vector{Float64}
         So_Max         :: Vector{Float64}
         θsMacMat_Max   ::	Vector{Float64}
         Φ_Max          ::	Vector{Float64}
		end # struct KOSUGI

		
		mutable struct VANGENUCHTEN # <>=<>=<>=<>=<>=<>=<>=<>=<>
         θs      ::	Vector{Float64}
         θr      ::	Vector{Float64}
         N       ::	Vector{Float64}
         Ψvg     ::	Vector{Float64}
         Ks      ::	Vector{Float64}
         Km      ::	Vector{Float64}
         Φ       ::	Vector{Float64}
			
         θs_Min  ::	Vector{Float64}
         θr_Min  ::	Vector{Float64}
         N_Min   ::	Vector{Float64}
         Ψvg_Min ::	Vector{Float64}
         Ks_Min  ::	Vector{Float64}
         Km_Min  ::	Vector{Float64}
         Φ_Min   ::	Vector{Float64}
			
         θs_Max  ::	Vector{Float64}
         θr_Max  ::	Vector{Float64}
         N_Max   ::	Vector{Float64}
         Ψvg_Max ::	Vector{Float64}
         Ks_Max  ::	Vector{Float64}
         Km_Max  ::	Vector{Float64}
         Φ_Max   ::	Vector{Float64}
		end # struct VANGENUCHTEN


		mutable struct BROOKS_COREY # <>=<>=<>=<>=<>=<>=<>=<>=<>
         θs      ::	Vector{Float64}
         θr      ::	Vector{Float64}
         λbc     ::	Vector{Float64}
         Ψbc     ::	Vector{Float64}
         Ks      ::	Vector{Float64}
         Φ       ::	Vector{Float64}
         Ψga     ::	Vector{Float64}
			
         θs_Min  ::	Vector{Float64}
         θr_Min  ::	Vector{Float64}
         λbc_Min ::	Vector{Float64}
         Ψbc_Min ::	Vector{Float64}
         Ks_Min  ::	Vector{Float64}
         Φ_Min   ::	Vector{Float64}
			
         θs_Max  ::	Vector{Float64}
         θr_Max  ::	Vector{Float64}
         λbc_Max ::	Vector{Float64}
         Ψbc_Max ::	Vector{Float64}
         Ks_Max  ::	Vector{Float64}
         Φ_Max   ::	Vector{Float64}
		end # struct BROOKS COREY


		mutable struct CLAPP_HORNBERGER # <>=<>=<>=<>=<>=<>=<>=<>=<>
			θs        ::	Vector{Float64}
         θr        ::	Vector{Float64}
         λch       ::	Vector{Float64}
         Ψch       ::	Vector{Float64}
         Ks        ::	Vector{Float64}
			Φ         ::	Vector{Float64}
         Ψga       ::	Vector{Float64}
			
			θs_Min        ::	Vector{Float64}
         θr_Min        ::	Vector{Float64}
         λch_Min       ::	Vector{Float64}
         Ψch_Min       ::	Vector{Float64}
         Ks_Min        ::	Vector{Float64}
			Φ_Min         ::	Vector{Float64}
			
			θs_Max        ::	Vector{Float64}
         θr_Max        ::	Vector{Float64}
         λch_Max       ::	Vector{Float64}
         Ψch_Max       ::	Vector{Float64}
         Ks_Max        ::	Vector{Float64}
         Φ_Max         ::	Vector{Float64}
		end # struct CLAPP_HORNBERGER


		mutable struct HYDRO_OTHER # <>=<>=<>=<>=<>=<>=<>=<>=<>
         Nse     :: 	Vector{Float64}
         Nse_θΨ  :: 	Vector{Float64}
         Nse_KΨ  :: 	Vector{Float64}
         Rmse    :: 	Vector{Float64}
         Rmse_θΨ :: 	Vector{Float64}
         Rmse_KΨ :: 	Vector{Float64}
		end

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROSTRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDROSTRUCT(optionₘ, N_iZ::Int64)
			# For all models
         θs         = fill(0.0::Float64, N_iZ)
         θr         = fill(0.0::Float64, N_iZ)
         Ks         = fill(0.0::Float64, N_iZ)
         θsMacMat   = fill(0.0::Float64, N_iZ)
         θsMacMat_ƞ = fill(0.0::Float64, N_iZ)
         Φ          = fill(0.0::Float64, N_iZ)
         So         = fill(0.0::Float64, N_iZ)

         θs_Min     = fill(0.0::Float64, N_iZ)
         θr_Min     = fill(0.0::Float64, N_iZ)
			Ks_Min     = fill(0.0::Float64, N_iZ)
			Φ_Min      = fill(0.0::Float64, N_iZ)
			So_Min     = fill(0.0::Float64, N_iZ)
			
         θs_Max     = fill(0.0::Float64, N_iZ)
         θr_Max     = fill(0.0::Float64, N_iZ)
			Ks_Max     = fill(0.0::Float64, N_iZ)
			Φ_Max      = fill(0.0::Float64, N_iZ)
			So_Max     = fill(0.0::Float64, N_iZ)
			
			if optionₘ.HydroModel⍰ == :Kosugi # <>=<>=<>=<>=<>
				σ              = fill(0.0::Float64, N_iZ)
            Ψm             = fill(0.0::Float64, N_iZ)
            σMac           = fill(0.0::Float64, N_iZ)
            ΨmMac          = fill(0.0::Float64, N_iZ)
				
            Ψm_Min         = fill(0.0::Float64, N_iZ)
            σ_Min          = fill(0.0::Float64, N_iZ)
            θsMacMat_ƞ_Min = fill(0.0::Float64, N_iZ)
            σMac_Min       = fill(0.0::Float64, N_iZ)
            ΨmMac_Min      = fill(0.0::Float64, N_iZ)
            θsMacMat_Min   = fill(0.0::Float64, N_iZ)

            Ψm_Max         = fill(0.0::Float64, N_iZ)
            σ_Max          = fill(0.0::Float64, N_iZ)
            θsMacMat_ƞ_Max = fill(0.0::Float64, N_iZ)
            σMac_Max       = fill(0.0::Float64, N_iZ)
            ΨmMac_Max      = fill(0.0::Float64, N_iZ)
            θsMacMat_Max   = fill(0.0::Float64, N_iZ)
          
				hydro = KOSUGI(θs, θr, Ks, σ, Ψm, θsMacMat_ƞ, σMac, ΨmMac, So, θsMacMat, Φ, θs_Min, θr_Min, Ks_Min, σ_Min,Ψm_Min, θsMacMat_ƞ_Min, σMac_Min, ΨmMac_Min, So_Min, θsMacMat_Min, Φ_Min, θs_Max, θr_Max, Ks_Max, σ_Max, Ψm_Max, θsMacMat_ƞ_Max, σMac_Max, ΨmMac_Max, So_Max, θsMacMat_Max, Φ_Max)
				return hydro

			elseif optionₘ.HydroModel⍰ == :Vangenuchten || optionₘ.HydroModel⍰ == :VangenuchtenJules # <>=<>=<>=<>=<>
            N       = fill(0.0::Float64, N_iZ)
            Ψvg     = fill(0.0::Float64, N_iZ)
            Km      = fill(0.0::Float64, N_iZ)

            θs_Min  = fill(0.0::Float64, N_iZ)
            θr_Min  = fill(0.0::Float64, N_iZ)
            N_Min   = fill(0.0::Float64, N_iZ)
            Ψvg_Min = fill(0.0::Float64, N_iZ)
            Ks_Min  = fill(0.0::Float64, N_iZ)
            Km_Min  = fill(0.0::Float64, N_iZ)
			
            θs_Max  = fill(0.0::Float64, N_iZ)
            θr_Max  = fill(0.0::Float64, N_iZ)
            N_Max   = fill(0.0::Float64, N_iZ)
            Ψvg_Max = fill(0.0::Float64, N_iZ)
            Ks_Max  = fill(0.0::Float64, N_iZ)
            Km_Max  = fill(0.0::Float64, N_iZ)

				hydro = VANGENUCHTEN(θs, θr, N, Ψvg, Ks, Km, Φ, θs_Min,θr_Min, N_Min, Ψvg_Min, Ks_Min, Km_Min, Φ_Min, θs_Max, θr_Max, N_Max, Ψvg_Max, Ks_Max, Km_Max, Φ_Max) 
				return hydro


			elseif optionₘ.HydroModel⍰ == :BrooksCorey # <>=<>=<>=<>=<>=<>
            λbc     = fill(0.0::Float64, N_iZ)
            Ψbc     = fill(0.0::Float64, N_iZ)
            λbc_Min = fill(0.0::Float64, N_iZ)
            Ψbc_Min = fill(0.0::Float64, N_iZ)
            λbc_Max = fill(0.0::Float64, N_iZ)
            Ψbc_Max = fill(0.0::Float64, N_iZ)
            Ψga     = fill(0.0::Float64, N_iZ)

				hydro = BROOKS_COREY(θs, θr, λbc, Ψbc, Ks, Φ, Ψga, θs_Min, θr_Min, λbc_Min, Ψbc_Min, Ks_Min, Φ_Min, θs_Max, θr_Max, λbc_Max, Ψbc_Max, Ks_Max, Φ_Max)
				return hydro

			elseif optionₘ.HydroModel⍰ == :ClappHornberger # <>=<>=<>=<>=<>=<>
				λch = fill(0.0::Float64, N_iZ)
				Ψch  = fill(0.0::Float64, N_iZ)
            Ψga  = fill(0.0::Float64, N_iZ)

				λch_Min = fill(0.0::Float64, N_iZ)
				Ψch_Min  = fill(0.0::Float64, N_iZ)

				λch_Max = fill(0.0::Float64, N_iZ)
				Ψch_Max  = fill(0.0::Float64, N_iZ)

				hydro = CLAPP_HORNBERGER(θs, θr, λch, Ψch, Ks, Φ, Ψga,  θs_Min, θr_Min, λch_Min, Ψch_Min, Ks_Min, Φ_Max, θs_Max, θr_Max, λch_Max, Ψch_Max, Ks_Max, Φ_Max)
				return hydro
			end # optionₘ.HydroModel⍰

		end #  function HYDROSTRUCT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYDRO_OTHER
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYDRO_OTHERS(N_iZ::Int64)
				Nse = fill(0.0::Float64, N_iZ)
            Nse_θΨ  = fill(0.0::Float64, N_iZ)
            Nse_KΨ  = fill(0.0::Float64, N_iZ)
            Rmse    = fill(0.0::Float64, N_iZ)
            Rmse_θΨ = fill(0.0::Float64, N_iZ)
            Rmse_KΨ = fill(0.0::Float64, N_iZ)

				hydroOther = HYDRO_OTHER(Nse, Nse_θΨ, Nse_KΨ, Rmse, Rmse_θΨ, Rmse_KΨ)	
				return hydroOther
			end  # function: HYDRO_OTHER
end  # module: hydroStruct

# ............................................................