# =============================================================
#		MODULE: name
# =============================================================
module hydroParam
	import ..option, ..param, ..psdThetar, ..hydroStruct
	export START_HYDROPARAM

		mutable struct OPTIMIZE
			Opt_θs :: Bool
			Opt_θr :: Bool
			Opt_Ks :: Bool
		end # struct VANGENUCHTEN

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_HYDROPARAM()
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_HYDROPARAM(N_SoilSelect, ∑Psd, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ)
			# INITIALIZING
				θ_Max 		= Array{Float64}(undef, (N_SoilSelect))
				θ_Min 		= Array{Float64}(undef, (N_SoilSelect))
				θr_Max 		= Array{Float64}(undef, (N_SoilSelect))
				θs_Min 		= Array{Float64}(undef, (N_SoilSelect))
				θs_Max 		= Array{Float64}(undef, (N_SoilSelect))
				K_KΨ_Max 	= Array{Float64}(undef, (N_SoilSelect))
				Ks_Min 		= Array{Float64}(undef, (N_SoilSelect))

			# WHAT TO OPTIMIZE (just initializing)
				Opt_θs 		= true
				Opt_θr 		= true
				Opt_Ks		= true

				opt = OPTIMIZE(Opt_θs, Opt_θr, Opt_Ks)

			# INITIALIZES HYDRAULIC PARAMETERS STRUCT INDEPENDENTLY OF THE SELECTED MODEL
				hydro = hydroStruct.HYDROSTRUCT(N_SoilSelect)

			# LOOPING FOR ERVERY SOIL
			for iSoil=1:N_SoilSelect
				# LIMITS
					θ_Max[iSoil] = maximum(θ_θΨ[iSoil, 1:N_θΨ[iSoil]])  	# Greatest measure θ
					θ_Min[iSoil] = minimum(θ_θΨ[iSoil, 1:N_θΨ[iSoil]])  	# Smallest measure θ

					if option.KunsatΨ
						K_KΨ_Max[iSoil] = maximum(K_KΨ[iSoil, 1:N_KΨ[iSoil]]) 	# Greatest measure of Kunsat
					end

				# DERIVING θr FROM DATA IF REQUESTED
					# Derive θr frpm PSD
					if (option.hydro.θrOpt == "Psd") && (option.Psd)
						hydro.θr[iSoil] = min( psdThetar.PSD_2_θr_FUNC(iSoil, ∑Psd), θ_Min[iSoil]-0.005 )
						opt.Opt_θr = false # No need to optimize θr
					# Keep θr = Cst
					elseif option.hydro.θrOpt == "Cst"
						hydro.θr[iSoil] = param.hydro.θr
						opt.Opt_θr = false
					# If optimised than maximum value of θr
					elseif option.hydro.θrOpt == "Opt"
						θr_Max[iSoil] = max( min(θ_Min[iSoil]-0.005, param.hydro.θr_Max), 0.0 ) # Maximum value of θr
						opt.Opt_θr = true
					end # option.hydro.θrOpt == "Psd"


				# DERIVING θs FROM DATA IF REQUESTED
					if option.hydro.θsOpt == "Data"
						hydro.θs[iSoil] = θ_Max[iSoil]
						opt.Opt_θs = false # No need to optimize θs
					elseif option.hydro.θsOpt == "Φ"
						hydro.θs[iSoil] = θ_Max[iSoil] * param.hydro.Coeff_Φ_2_θs
						opt.Opt_θs = false
					elseif option.hydro.θsOpt == "Opt"
						θs_Min[iSoil] = θ_Max[iSoil]
						θs_Max[iSoil] = θ_Max[iSoil] * param.hydro.Coeff_θs_Max
						opt.Opt_θs = true
					end # option.hydro.θsOpt


				# DERIVING Ks FROM DATA IF REQUESTED
					if option.hydro.KsOpt == "Data" && option.KunsatΨ
						hydro.Ks[iSoil] = K_KΨ_Max[iSoil]
						opt.Opt_Ks = false
					elseif option.hydro.KsOpt == "Opt" && option.KunsatΨ
						Ks_Min[iSoil] = K_KΨ_Max[iSoil]
						opt.Opt_Ks = true
					end # option.KunsatΨ
			end  # for iSoil=1:N_SoilSelect

			# DETERMENING THE NUMBER OF PARAMETERS TO BE OPTIMIZED
				if option.hydro.UnimodalBimodal == "Unimodal"
					N_ParamOpt = 5 	# Number of parameters to be optimized (will change)
				elseif option.hydro.UnimodalBimodal == "Bimodal"
					N_ParamOpt = 8 	# Number of parameters to be optimized (will change)
				end  # if: option.hydro.UnimodalBimodal
				if !opt.Opt_θr
					N_ParamOpt -= 1
				end
				if !opt.Opt_θs
					 N_ParamOpt -= 1
				end
				if !opt.Opt_Ks
					N_ParamOpt -= 1
				end 
				println("	...Optimizing $N_ParamOpt hydraulic parameters")  
				
			
			# OPTIMIZATION HYDRAULIC PARAMETERS
				if option.hydro.HydroModel == "Kosugi"
					Of, Of_θΨ, Of_Kunsat, hydro = kg.HYDROPARAM_OPT(θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, θr_Max, θs_Min, θs_Max, Ks_Min, N_ParamOpt, N_SoilSelect, hydro, opt)
					
				elseif option.hydro.HydroModel == "Vangenuchten"
					Of, Of_θΨ, Of_Kunsat, hydro = vg.HYDROPARAM_OPT(θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, θr_Max, θs_Min, θs_Max, Ks_Min, N_ParamOpt, N_SoilSelect, hydro, opt)
				end

			return Of, Of_θΨ, Of_Kunsat, hydro
		end # function: START_HYDROPARAM

		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

		# =============================================================
		#		MODULE: kg
		# =============================================================
		module kg
			import ...ofHydro, ...param, ..option
			using BlackBoxOptim
			export HYDROPARAM_OPT
	
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : HYDROPARAM_OPT
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function HYDROPARAM_OPT(θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, θr_Max, θs_Min, θs_Max, Ks_Min, N_ParamOpt, N_SoilSelect, hydro, opt)
					Of 			= Array{Float64}(undef, (N_SoilSelect))
					Of_θΨ 		= Array{Float64}(undef, (N_SoilSelect))
					Of_Kunsat 	= zeros(Float64, (N_SoilSelect))

					for iSoil=1:N_SoilSelect

						if option.hydro.UnimodalBimodal=="Bimodal" && opt.Opt_θs && opt.Opt_θr && opt.Opt_Ks

							SearchRanges =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (0.0, θr_Max[iSoil]), (θs_Min[iSoil], θs_Max[iSoil]), (param.hydro.∇_θsMat_Min, 1.0), (0.0, 1.0), (log10(param.hydro.ΨmMac_Min), log10(param.hydro.ΨmMac_Max)), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θr=P[3], θs=P[4], ∇_θsMat=P[5], ∇_σMac=P[6], ΨmMac=10.0^P[7], Ks=P[8])[1]; SearchRange = SearchRanges, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[4]
							∇_θsMat 			= BlackBoxOptim.best_candidate(Optimization)[5]
							∇_σMac 				= BlackBoxOptim.best_candidate(Optimization)[6]
							hydro.ΨmMac[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[7]
							hydro.Ks[iSoil] 	=  BlackBoxOptim.best_candidate(Optimization)[8]

							hydro.θsMat[iSoil] = ∇NORM_2_PARAMETER(∇_θsMat, hydro.θr[iSoil], hydro.θs[iSoil])
							hydro.σMac[iSoil]  = ∇NORM_2_PARAMETER(∇_σMac, param.hydro.σMac_Min, hydro.σ[iSoil])


						elseif option.hydro.UnimodalBimodal=="Bimodal" && opt.Opt_θs && !opt.Opt_θr && opt.Opt_Ks

							SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (θs_Min[iSoil], θs_Max[iSoil]), (param.hydro.∇_θsMat_Min, 1.0), (0.0, 1.0), (log10(param.hydro.ΨmMac_Min), log10(param.hydro.ΨmMac_Max)), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θs=P[3], ∇_θsMat=P[4], ∇_σMac=P[5], ΨmMac=10.0^P[6], Ks=P[7])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							∇_θsMat 			= BlackBoxOptim.best_candidate(Optimization)[4]
							∇_σMac 				= BlackBoxOptim.best_candidate(Optimization)[5]
							hydro.ΨmMac[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[6]
							hydro.Ks[iSoil] 	=  BlackBoxOptim.best_candidate(Optimization)[7]
							hydro.θsMat[iSoil] = ∇NORM_2_PARAMETER(∇_θsMat, hydro.θr[iSoil], hydro.θs[iSoil])
							hydro.σMac[iSoil]  = ∇NORM_2_PARAMETER(∇_σMac, param.hydro.σMac_Min, hydro.σ[iSoil])


						elseif option.hydro.UnimodalBimodal=="Bimodal" && !opt.Opt_θs && opt.Opt_θr && opt.Opt_Ks

							SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (0.0, θr_Max[iSoil]), (param.hydro.∇_θsMat_Min, 1.0), (0.0, 1.0), (log10(param.hydro.ΨmMac_Min), log10(param.hydro.ΨmMac_Max)), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θr=P[3], ∇_θsMat=P[4], ∇_σMac=P[5], ΨmMac=10.0^P[6], Ks=P[7])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							∇_θsMat 			= BlackBoxOptim.best_candidate(Optimization)[4]
							∇_σMac 				= BlackBoxOptim.best_candidate(Optimization)[5]
							hydro.ΨmMac[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[6]
							hydro.Ks[iSoil] 	=  BlackBoxOptim.best_candidate(Optimization)[7]
							hydro.θsMat[iSoil] = ∇NORM_2_PARAMETER(∇_θsMat, hydro.θr[iSoil], hydro.θs[iSoil])
							hydro.σMac[iSoil]  = ∇NORM_2_PARAMETER(∇_σMac, param.hydro.σMac_Min, hydro.σ[iSoil])


						elseif option.hydro.UnimodalBimodal=="Bimodal" && !opt.Opt_θs && !opt.Opt_θr && opt.Opt_Ks

							SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (param.hydro.∇_θsMat_Min, 1.0), (0.0, 1.0), (log10(param.hydro.ΨmMac_Min), log10(param.hydro.ΨmMac_Max)), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], ∇_θsMat=P[3], ∇_σMac=P[4], ΨmMac=10.0^P[5], Ks=P[6])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							∇_θsMat 			= BlackBoxOptim.best_candidate(Optimization)[3]
							∇_σMac 				= BlackBoxOptim.best_candidate(Optimization)[4]
							hydro.ΨmMac[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[5]
							hydro.Ks[iSoil] 	=  BlackBoxOptim.best_candidate(Optimization)[6]
							hydro.θsMat[iSoil] = ∇NORM_2_PARAMETER(∇_θsMat, hydro.θr[iSoil], hydro.θs[iSoil])
							hydro.σMac[iSoil]  = ∇NORM_2_PARAMETER(∇_σMac, param.hydro.σMac_Min, hydro.σ[iSoil])


						elseif option.hydro.UnimodalBimodal=="Unimodal" && opt.Opt_θs && opt.Opt_θr && opt.Opt_Ks

							SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (0.0, θr_Max[iSoil]), (θs_Min[iSoil], θs_Max[iSoil]), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θr=P[3], θs=P[4], Ks=P[5])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[4]
							hydro.Ks[iSoil] 	=  BlackBoxOptim.best_candidate(Optimization)[5]
							hydro.θsMat[iSoil] = hydro.θs[iSoil]
							hydro.σMac[iSoil]  = hydro.σ[iSoil]


						elseif option.hydro.UnimodalBimodal=="Unimodal" && opt.Opt_θs && !opt.Opt_θr && opt.Opt_Ks

							SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (θs_Min[iSoil], θs_Max[iSoil]), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θs=P[3], Ks=P[4])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.Ks[iSoil] 	=  BlackBoxOptim.best_candidate(Optimization)[4]
							hydro.θsMat[iSoil] = hydro.θs[iSoil]
							hydro.σMac[iSoil]  = hydro.σ[iSoil]


						elseif option.hydro.UnimodalBimodal=="Unimodal" && !opt.Opt_θs && opt.Opt_θr && opt.Opt_Ks
							#############
							SearchRanges =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (0.0, θr_Max[iSoil]), (log10(Ks_Min[iSoil]), log10(param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θr=P[3], Ks=10.0^P[4])[1]; SearchRange = SearchRanges, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.Ks[iSoil] 	=  10.0 ^ BlackBoxOptim.best_candidate(Optimization)[4]
							hydro.θsMat[iSoil] = hydro.θs[iSoil]
							hydro.σMac[iSoil]  = hydro.σ[iSoil]


						elseif option.hydro.UnimodalBimodal=="Unimodal" && !opt.Opt_θs && !opt.Opt_θr && opt.Opt_Ks

							SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], Ks=P[3])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.Ks[iSoil] 	=  BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.θsMat[iSoil] = hydro.θs[iSoil]
							hydro.σMac[iSoil]  = hydro.σ[iSoil]


							elseif option.hydro.UnimodalBimodal=="Bimodal" && opt.Opt_θs && opt.Opt_θr && !opt.Opt_Ks

								SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (0.0, θr_Max[iSoil]), (θs_Min[iSoil], θs_Max[iSoil]), (param.hydro.∇_θsMat_Min, 1.0), (0.0, 1.0), (log10(param.hydro.ΨmMac_Min), log10(param.hydro.ΨmMac_Max))]
	
								Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θr=P[3], θs=P[4], ∇_θsMat=P[5], ∇_σMac=P[6], ΨmMac=10.0^P[7])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
			
								hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
								hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
								hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
								hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[4]
								∇_θsMat 			= BlackBoxOptim.best_candidate(Optimization)[5]
								∇_σMac 				= BlackBoxOptim.best_candidate(Optimization)[6]
								hydro.ΨmMac[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[7]
								hydro.θsMat[iSoil] = ∇NORM_2_PARAMETER(∇_θsMat, hydro.θr[iSoil], hydro.θs[iSoil])
								hydro.σMac[iSoil]  = ∇NORM_2_PARAMETER(∇_σMac, param.hydro.σMac_Min, hydro.σ[iSoil])
		
		
							elseif option.hydro.UnimodalBimodal=="Bimodal" && opt.Opt_θs && !opt.Opt_θr && !opt.Opt_Ks
	
								SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (θs_Min[iSoil], θs_Max[iSoil]), (param.hydro.∇_θsMat_Min, 1.0), (0.0, 1.0), (log10(param.hydro.ΨmMac_Min), log10(param.hydro.ΨmMac_Max))]
	
								Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θs=P[3], ∇_θsMat=P[4], ∇_σMac=P[5], ΨmMac=10.0^P[6])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
			
								hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
								hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
								hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
								∇_θsMat 			= BlackBoxOptim.best_candidate(Optimization)[4]
								∇_σMac 				= BlackBoxOptim.best_candidate(Optimization)[5]
								hydro.ΨmMac[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[6]
								hydro.θsMat[iSoil] = ∇NORM_2_PARAMETER(∇_θsMat, hydro.θr[iSoil], hydro.θs[iSoil])
								hydro.σMac[iSoil]  = ∇NORM_2_PARAMETER(∇_σMac, param.hydro.σMac_Min, hydro.σ[iSoil])
	
	
							elseif option.hydro.UnimodalBimodal=="Bimodal" && !opt.Opt_θs && opt.Opt_θr && !opt.Opt_Ks
	
								SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (0.0, θr_Max[iSoil]), (param.hydro.∇_θsMat_Min, 1.0), (0.0, 1.0), (log10(param.hydro.ΨmMac_Min), log10(param.hydro.ΨmMac_Max))]
	
								Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θr=P[3], ∇_θsMat=P[4], ∇_σMac=P[5], ΨmMac=10.0^P[6])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
			
								hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
								hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
								hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
								∇_θsMat 			= BlackBoxOptim.best_candidate(Optimization)[4]
								∇_σMac 				= BlackBoxOptim.best_candidate(Optimization)[5]
								hydro.ΨmMac[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[6]
								hydro.θsMat[iSoil] = ∇NORM_2_PARAMETER(∇_θsMat, hydro.θr[iSoil], hydro.θs[iSoil])
								hydro.σMac[iSoil]  = ∇NORM_2_PARAMETER(∇_σMac, param.hydro.σMac_Min, hydro.σ[iSoil])
	
	
							elseif option.hydro.UnimodalBimodal=="Bimodal" && !opt.Opt_θs && !opt.Opt_θr && !opt.Opt_Ks
	
								SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (param.hydro.∇_θsMat_Min, 1.0), (0.0, 1.0), (log10(param.hydro.ΨmMac_Min), log10(param.hydro.ΨmMac_Max)), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]
	
								Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], ∇_θsMat=P[3], ∇_σMac=P[4], ΨmMac=10.0^P[5])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
			
								hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
								hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
								∇_θsMat 			= BlackBoxOptim.best_candidate(Optimization)[3]
								∇_σMac 				= BlackBoxOptim.best_candidate(Optimization)[4]
								hydro.ΨmMac[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[5]
								hydro.θsMat[iSoil] = ∇NORM_2_PARAMETER(∇_θsMat, hydro.θr[iSoil], hydro.θs[iSoil])
								hydro.σMac[iSoil]  = ∇NORM_2_PARAMETER(∇_σMac, param.hydro.σMac_Min, hydro.σ[iSoil])
	
	
							elseif option.hydro.UnimodalBimodal=="Unimodal" && opt.Opt_θs && opt.Opt_θr && !opt.Opt_Ks
	
									SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (0.0, θr_Max[iSoil]), (θs_Min[iSoil], θs_Max[iSoil])]
		
									Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θr=P[3], θs=P[4])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
				
									hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
									hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
									hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
									hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[4]
									hydro.θsMat[iSoil] = hydro.θs[iSoil]
									hydro.σMac[iSoil]  = hydro.σ[iSoil]
	
	
								elseif option.hydro.UnimodalBimodal=="Unimodal" && opt.Opt_θs && !opt.Opt_θr && !opt.Opt_Ks
	
									SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (θs_Min[iSoil], θs_Max[iSoil])]
		
									Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θs=P[3])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
				
									hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
									hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
									hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
									hydro.θsMat[iSoil] = hydro.θs[iSoil]
									hydro.σMac[iSoil]  = hydro.σ[iSoil]
	
	
								elseif option.hydro.UnimodalBimodal=="Unimodal" && !opt.Opt_θs && opt.Opt_θr && !opt.Opt_Ks
	
									SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max)), (0.0, θr_Max[iSoil])]
		
									Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θr=P[3])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
				
									hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
									hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
									hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
									hydro.θsMat[iSoil] = hydro.θs[iSoil]
									hydro.σMac[iSoil]  = hydro.σ[iSoil]
	
	
								elseif option.hydro.UnimodalBimodal=="Unimodal" && !opt.Opt_θs && !opt.Opt_θr && !opt.Opt_Ks
	
									SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (log10(param.hydro.Ψm_Min), log10(param.hydro.Ψm_Max))]
		
									Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
				
									hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
									hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
									hydro.θsMat[iSoil] = hydro.θs[iSoil]
									hydro.σMac[iSoil]  = hydro.σ[iSoil]

							end #Option

	
						Of[iSoil], Of_θΨ[iSoil], Of_Kunsat[iSoil] = OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro)
	
					end  # for iSoil=1:N_SoilSelect
					
					return Of, Of_θΨ, Of_Kunsat, hydro
				end  # function: HYDROPARAM_OPT


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : HYDRO_PROCESS
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=hydro.σ[iSoil], Ψm=hydro.Ψm[iSoil], θr=hydro.θr[iSoil], θs=hydro.θs[iSoil], Ks=hydro.Ks[iSoil], ∇_θsMat=1.0, ∇_σMac=1.0, ΨmMac=hydro.Ψm[iSoil])
	
					hydro.θsMat[iSoil] = ∇NORM_2_PARAMETER(∇_θsMat, θr, θs)
					
					hydro.σMac[iSoil]  = ∇NORM_2_PARAMETER(∇_σMac, param.hydro.σMac_Min, σ)

					hydro.θs[iSoil] = θs
					hydro.θr[iSoil] = θr
					hydro.Ks[iSoil] = Ks
					hydro.σ[iSoil] = σ
					hydro.Ψm[iSoil] = Ψm
					hydro.ΨmMac[iSoil] = ΨmMac

					Of, Of_θΨ, Of_Kunsat = ofHydro.OF_WRC_KUNSAT(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro) 
	
					return Of, Of_θΨ, Of_Kunsat
				end  # function: HYDRO_PROCESS
					

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : HYDRO_ADJEUSTMENTS
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function ∇NORM_2_PARAMETER(∇P, P_Min, P_Max)
					return P = ∇P * (P_Max - P_Min) + P_Min
				end  # function: HYDRO_ADJEUSTMENTS
			
		end  # module kg
		# ............................................................


		# =============================================================
		#		MODULE: vg
		# =============================================================
		module vg
		import ...ofHydro, ...param, ..option
		using BlackBoxOptim
		export HYDROPARAM_OPT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		vg FUNCTION : HYDROPARAM_OPT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYDROPARAM_OPT(θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, θr_Max, θs_Min, θs_Max, Ks_Min, N_ParamOpt, N_SoilSelect, hydro, opt)
				Of 			= Array{Float64}(undef, (N_SoilSelect))
				Of_θΨ 		= Array{Float64}(undef, (N_SoilSelect))
				Of_Kunsat 	= Array{Float64}(undef, (N_SoilSelect))

				for iSoil=1:N_SoilSelect

					# it is to be noted that option Opt_Ks is not needed since it is regulated by N_ParamOpt
					if option.hydro.UnimodalBimodal=="Unimodal" && opt.Opt_θs && opt.Opt_θr

							SearchRange =[(param.hydro.N_Min, param.hydro.N_Max), (log10(param.hydro.Ψvg_Min), log10(param.hydro.Ψvg_Max)), (0.0, θr_Max[iSoil]), (θs_Min[iSoil], θs_Max[iSoil]), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; N=P[1], Ψvg=10.0^P[2], θr=P[3], θs=P[4], Ks=P[5])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.N[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψvg[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[4]
							hydro.Ks[iSoil] 	=  BlackBoxOptim.best_candidate(Optimization)[5]


						elseif option.hydro.UnimodalBimodal=="Unimodal" && opt.Opt_θs && !opt.Opt_θr

							SearchRange =[(param.hydro.N_Min, param.hydro.N_Max), (log10(param.hydro.Ψvg_Min), log10(param.hydro.Ψvg_Max)), (θs_Min[iSoil], θs_Max[iSoil]), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; N=P[1], Ψvg=10.0^P[2], θs=P[3], Ks=P[4])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.N[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψvg[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.Ks[iSoil] 	=  BlackBoxOptim.best_candidate(Optimization)[4]


						elseif option.hydro.UnimodalBimodal=="Unimodal" && !opt.Opt_θs && opt.Opt_θr

							SearchRange =[(param.hydro.N_Min, param.hydro.N_Max), (log10(param.hydro.Ψvg_Min), log10(param.hydro.Ψvg_Max)), (0.0, θr_Max[iSoil]), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; N=P[1], Ψvg=10.0^P[2], θr=P[3], Ks=P[4])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.N[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψvg[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.Ks[iSoil] 	=  BlackBoxOptim.best_candidate(Optimization)[4]


						elseif option.hydro.UnimodalBimodal=="Unimodal" && !opt.Opt_θs && !opt.Opt_θr

							SearchRange =[(param.hydro.N_Min, param.hydro.N_Max), (log10(param.hydro.Ψvg_Min), log10(param.hydro.Ψvg_Max)), ((Ks_Min[iSoil]), (param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; N=P[1], Ψvg=10.0^P[2], Ks=P[3])[1]; SearchRange = SearchRange, NumDimensions=N_ParamOpt, TraceMode=:silent)
		
							hydro.N[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψvg[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.Ks[iSoil] 	=  BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.θsMat[iSoil] = hydro.θs[iSoil]
							hydro.σMac[iSoil]  = hydro.N[iSoil]

						end #Option

					Of[iSoil], Of_θΨ[iSoil], Of_Kunsat[iSoil] = OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro)

				end  # for iSoil=1:N_SoilSelect
				
				return Of, Of_θΨ, Of_Kunsat, hydro
			end  # function: HYDROPARAM_OPT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		vg FUNCTION : HYDRO_PROCESS
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; N=hydro.N[iSoil], Ψvg=hydro.Ψvg[iSoil], θr=hydro.θr[iSoil], θs=hydro.θs[iSoil], Ks=hydro.Ks[iSoil])

				hydro.θs[iSoil] = θs
				hydro.θr[iSoil] = θr
				hydro.Ks[iSoil] = Ks
				hydro.N[iSoil] = N
				hydro.Ψvg[iSoil] = Ψvg

				Of, Of_θΨ, Of_Kunsat = ofHydro.OF_WRC_KUNSAT(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro) 

				return Of, Of_θΨ, Of_Kunsat
			end  # function: HYDRO_PROCESS
			
		end  # module vg
		# ............................................................
		
end  # module: hydro