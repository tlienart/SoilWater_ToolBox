# =============================================================
#		MODULE: module hydroInitialize

# =============================================================
module hydroInitialize
	import ..option, ..param, ..psdThetar, ..hydroStruct
	export HYDRO_INITIALIZE

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDRO_INITIALIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	mutable struct OPTIMIZE
        Opt_θs     :: Bool
        Opt_θr     :: Bool
        Opt_Ks     :: Bool
        N_ParamOpt :: Int
	end # struct 

	function HYDRO_INITIALIZE(N_SoilSelect, ∑Psd, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro, optionHydro)

		# INITIALIZING
			θ_Max 		=  Array{Float64}(undef, (N_SoilSelect))
			θ_Min 		=  Array{Float64}(undef, (N_SoilSelect))
			θr_Max 		=  Array{Float64}(undef, (N_SoilSelect)) 
			θs_Min 		=  Array{Float64}(undef, (N_SoilSelect))
			θs_Max 		= Array{Float64}(undef, (N_SoilSelect))
			K_KΨ_Max 	= zeros(Float64, N_SoilSelect)
			Ks_Min 		= zeros(Float64, N_SoilSelect)

			Opt_θs 		= true
			Opt_θr 		= true
			Opt_Ks		= true
			N_ParamOpt	= 1

			opt = OPTIMIZE(Opt_θs, Opt_θr, Opt_Ks, N_ParamOpt)

		# LOOPING FOR ERVERY SOIL
			for iSoil=1:N_SoilSelect
				# LIMITS
					θ_Min[iSoil] = minimum(θ_θΨ[iSoil, 1:N_θΨ[iSoil]])  	# Smallest measure θ
					θ_Max[iSoil] = maximum(θ_θΨ[iSoil, 1:N_θΨ[iSoil]])  	# Greatest measure θ

					if optionHydro.KunsatΨ
						K_KΨ_Max[iSoil] = maximum(K_KΨ[iSoil, 1:N_KΨ[iSoil]]) 	# Greatest measure of Kunsat
					end

				# DERIVING θr FROM DATA IF REQUESTED
					
					if (optionHydro.θrOpt == "Psd") && (option.Psd) # Derive θr frpm PSD
						hydro.θr[iSoil] = min( psdThetar.PSD_2_θr_FUNC(iSoil, ∑Psd), θ_Min[iSoil]-0.005 )
						opt.Opt_θr = false # No need to optimize θr
				
					elseif optionHydro.θrOpt == "Cst" # Keep θr = Cst<>=<>=<>=<>=<>
						hydro.θr[iSoil] = param.hydro.θr
						opt.Opt_θr = false
			
					elseif optionHydro.θrOpt == "Opt" # <>=<>=<>=<>=<>
						θr_Max[iSoil] = max( min(θ_Min[iSoil]-0.005, param.hydro.θr_Max), 0.0 ) # Maximum value of θr
						opt.Opt_θr = true

					elseif optionHydro.θrOpt == "Known" # <>=<>=<>=<>=<>
						opt.Opt_θr = false

					end # optionHydro.θrOpt == "Psd"


				# DERIVING θs FROM DATA IF REQUESTED
					if optionHydro.θsOpt == "Data" # TODO need to derive from bulk density 
						hydro.θs[iSoil] = θ_Max[iSoil]
						opt.Opt_θs = false # No need to optimize θs
					
					elseif optionHydro.θsOpt == "Φ" # <>=<>=<>=<>=<>
						hydro.θs[iSoil] = max(θ_Max[iSoil] * param.hydro.Coeff_Φ_2_θs, θ_θΨ[iSoil, 2] + 0.005) # So that monotically increasing
						opt.Opt_θs = false

					elseif optionHydro.θsOpt == "Opt" # <>=<>=<>=<>=<>
						θs_Min[iSoil] = θ_θΨ[iSoil, 2] + 0.005
						θs_Max[iSoil] = max(θ_Max[iSoil] * param.hydro.Coeff_θs_Max, θ_θΨ[iSoil, 2] + 0.005)
						opt.Opt_θs = true

					elseif optionHydro.θsOpt == "Known" # <>=<>=<>=<>=<>
						opt.Opt_θs = false

					end # optionHydro.θsOpt


				# DERIVING Ks FROM DATA IF REQUESTED
					if !(optionHydro.KunsatΨ)
						opt.Opt_Ks = false

					elseif optionHydro.KsOpt == "Data" && optionHydro.KunsatΨ
						hydro.Ks[iSoil] = K_KΨ_Max[iSoil]
						opt.Opt_Ks = false

					elseif optionHydro.KsOpt == "Opt" && optionHydro.KunsatΨ  # <>=<>=<>=<>=<>
						Ks_Min[iSoil] = K_KΨ_Max[iSoil]
						opt.Opt_Ks = true

					elseif optionHydro.KsOpt == "Known" # <>=<>=<>=<>=<>
						opt.Opt_Ks = false

					end # optionHydro.KunsatΨ
			end  # for iSoil=1:N_SoilSelect

	# DETERMENING THE NUMBER OF PARAMETERS TO BE OPTIMIZED
		if optionHydro.UnimodalBimodal == "Unimodal"
			opt.N_ParamOpt = 5 	# Number of parameters to be optimized (will change)
		elseif optionHydro.UnimodalBimodal == "Bimodal"
			opt.N_ParamOpt = 8 	# Number of parameters to be optimized (will change)
		end  # if: optionHydro.UnimodalBimodal
		if !opt.Opt_θr
			opt.N_ParamOpt -= 1
		end
		if !opt.Opt_θs
			opt.N_ParamOpt -= 1
		end
		if !opt.Opt_Ks
			opt.N_ParamOpt -= 1
		end 
		println("    ~ Optimizing  $(opt.N_ParamOpt)  hydraulic parameters  ~")  
		
		return opt, θr_Max, θs_Min, θs_Max, Ks_Min, hydro
	end  # function: HYDRO_INITIALIZE
	
end  # module hydroInitialize

# ............................................................