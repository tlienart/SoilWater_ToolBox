# =============================================================
#		MODULE: name
# =============================================================
module hydroParam
	import ..option, ..param, ..hydroInitialize
	export START_HYDROPARAM

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_HYDROPARAM()
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_HYDROPARAM(N_SoilSelect, ∑Psd, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ)

			# INITIALIZATION
				opt, θr_Max, θs_Min, θs_Max, Ks_Min, hydro = hydroInitialize.HYDRO_INITIALIZE(N_SoilSelect, ∑Psd, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ)
			
			# OPTIMIZATION HYDRAULIC PARAMETERS
				if option.hydro.HydroModel == "Kosugi" # <>=<>=<>=<>=<>
					hydro = kg.HYDROPARAM_OPT(θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, θr_Max, θs_Min, θs_Max, Ks_Min, N_SoilSelect, hydro, opt)
					
				elseif option.hydro.HydroModel == "Vangenuchten" # <>=<>=<>=<>=<>
					hydro = vg.HYDROPARAM_OPT(θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, θr_Max, θs_Min, θs_Max, Ks_Min, N_SoilSelect, hydro, opt)
				end

			return hydro
		end # function: START_HYDROPARAM

		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

		# =============================================================
		#		MODULE: kg
		# =============================================================
		module kg
			import ...ofHydro, ...param, ..option
			import BlackBoxOptim
			export HYDROPARAM_OPT
	
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : HYDROPARAM_OPT
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function HYDROPARAM_OPT(θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, θr_Max, θs_Min, θs_Max, Ks_Min, N_SoilSelect, hydro, opt)

					for iSoil=1:N_SoilSelect
						if option.hydro.UnimodalBimodal=="Bimodal" && !opt.Opt_θs && opt.Opt_θr && opt.Opt_Ks  # This one <>=<>=<>=<>=<>

							SearchRange =[(param.hydro.kg.σ_Min, param.hydro.kg.σ_Max), (log10(param.hydro.kg.Ψm_Min), log10(param.hydro.kg.Ψm_Max)), (0.0, θr_Max[iSoil]), (param.hydro.kg.∇_θsMat_Min, param.hydro.kg.∇_θsMat_Max), (param.hydro.kg.σMac_Min, param.hydro.kg.σMac_Max), (param.hydro.kg.ΨmMac_Min, param.hydro.kg.ΨmMac_Max), (log10(Ks_Min[iSoil]), log10(param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θr=P[3], ∇_θsMat=P[4], σMac=P[5], ΨmMac=10.0^P[6], Ks=10.0^P[7])[1]; SearchRange = SearchRange, NumDimensions=opt.N_ParamOpt, TraceMode=:silent)
		
							hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							∇_θsMat 			= BlackBoxOptim.best_candidate(Optimization)[4]
							hydro.σMac[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[5]
							hydro.ΨmMac[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[6]
							hydro.Ks[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[7]

							hydro.θsMat[iSoil] = ∇NORM_2_PARAMETER(∇_θsMat, hydro.θr[iSoil], hydro.θs[iSoil])

						elseif option.hydro.UnimodalBimodal=="Unimodal" && opt.Opt_θs && opt.Opt_θr && opt.Opt_Ks  # <>=<>=<>=<>=<>

							SearchRange =[(param.hydro.kg.σ_Min, param.hydro.kg.σ_Max), (log10(param.hydro.kg.Ψm_Min), log10(param.hydro.kg.Ψm_Max)), (0.0, θr_Max[iSoil]), (θs_Min[iSoil], θs_Max[iSoil]), (log10(Ks_Min[iSoil]), log10(param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=P[1], Ψm=10.0^P[2], θr=P[3], θs=P[4], Ks=10.0^P[5])[1]; SearchRange = SearchRange, NumDimensions=opt.N_ParamOpt, TraceMode=:silent)
		
							hydro.σ[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψm[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[4]
							hydro.Ks[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[5]
							hydro.θsMat[iSoil] = hydro.θs[iSoil]
							hydro.σMac[iSoil]  = hydro.σ[iSoil]

							Of = BlackBoxOptim.best_fitness(Optimization)

						else
							error( " SoilWater-Toolbox error: option.hydro not found ")
						end #Option

						# STATISTICS OF THE OPTIMIZATION
							Of, Of_θΨ, Of_Kunsat = ofHydro.OF_WRC_KUNSAT(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro) 

							hydro.Nse[iSoil] = 1.0 - Of
							hydro.Nse_θψ[iSoil] = 1.0 - Of_θΨ
							hydro.Nse_Kψ[iSoil] = 1.0 - Of_Kunsat
	
					end  # for iSoil=1:N_SoilSelect
					
					return hydro
				end  # function: HYDROPARAM_OPT


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : HYDRO_PROCESS
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=hydro.σ[iSoil], Ψm=hydro.Ψm[iSoil], θr=hydro.θr[iSoil], θs=hydro.θs[iSoil], Ks=hydro.Ks[iSoil], ∇_θsMat=1.0, σMac=hydro.σMac[iSoil], ΨmMac=hydro.ΨmMac[iSoil])

					hydro.θs[iSoil] = θs
					hydro.θr[iSoil] = θr
					hydro.Ks[iSoil] = Ks
					hydro.σ[iSoil] = σ
					hydro.Ψm[iSoil] = Ψm
					hydro.ΨmMac[iSoil] = ΨmMac
					hydro.σMac[iSoil] = σMac
	
					hydro.θsMat[iSoil] = ∇NORM_2_PARAMETER(∇_θsMat, hydro.θr[iSoil], hydro.θs[iSoil])

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
		import BlackBoxOptim
		export HYDROPARAM_OPT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		vg FUNCTION : HYDROPARAM_OPT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYDROPARAM_OPT(θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, θr_Max, θs_Min, θs_Max, Ks_Min, N_SoilSelect, hydro, opt)

				for iSoil=1:N_SoilSelect

					# it is to be noted that option Opt_Ks is not needed since it is regulated by opt.N_ParamOpt
					if option.hydro.UnimodalBimodal=="Unimodal" && opt.Opt_θs && opt.Opt_θr  # <>=<>=<>=<>=<>

							SearchRange =[(param.hydro.vg.N_Min, param.hydro.vg.N_Max), (log10(param.hydro.vg.Ψvg_Min), log10(param.hydro.vg.Ψvg_Max)), (0.0, θr_Max[iSoil]), (θs_Min[iSoil], θs_Max[iSoil]), (log10(Ks_Min[iSoil]), log10(param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; N=P[1], Ψvg=10.0^P[2], θr=P[3], θs=P[4], Ks=10.0^P[5])[1]; SearchRange = SearchRange, NumDimensions=opt.N_ParamOpt, TraceMode=:silent)
		
							hydro.N[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψvg[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[4]
							hydro.Ks[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[5]


						elseif option.hydro.UnimodalBimodal=="Unimodal" && opt.Opt_θs && !opt.Opt_θr  # <>=<>=<>=<>=<>

							SearchRange =[(param.hydro.vg.N_Min, param.hydro.vg.N_Max), (log10(param.hydro.vg.Ψvg_Min), log10(param.hydro.vg.Ψvg_Max)), (θs_Min[iSoil], θs_Max[iSoil]), (log10(Ks_Min[iSoil]), log10(param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; N=P[1], Ψvg=10.0^P[2], θs=P[3], Ks=10.0^P[4])[1]; SearchRange = SearchRange, NumDimensions=opt.N_ParamOpt, TraceMode=:silent)
		
							hydro.N[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψvg[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θs[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.Ks[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[4]


						elseif option.hydro.UnimodalBimodal=="Unimodal" && !opt.Opt_θs && opt.Opt_θr  # test <>=<>=<>=<>=<>

							SearchRange =[(param.hydro.vg.N_Min, param.hydro.vg.N_Max), (log10(param.hydro.vg.Ψvg_Min), log10(param.hydro.vg.Ψvg_Max)), (0.0, θr_Max[iSoil]), (log10(Ks_Min[iSoil]), log10(param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; N=P[1], Ψvg=10.0^P[2], θr=P[3],Ks= 10.0^P[4])[1]; SearchRange = SearchRange, NumDimensions=opt.N_ParamOpt, TraceMode=:silent)
		
							hydro.N[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψvg[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.θr[iSoil] 	= BlackBoxOptim.best_candidate(Optimization)[3]
							hydro.Ks[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[4]


						elseif option.hydro.UnimodalBimodal=="Unimodal" && !opt.Opt_θs && !opt.Opt_θr  # <>=<>=<>=<>=<>

							SearchRange =[(param.hydro.vg.N_Min, param.hydro.vg.N_Max), (log10(param.hydro.vg.Ψvg_Min), log10(param.hydro.vg.Ψvg_Max)), (log10(Ks_Min[iSoil]), log10(param.hydro.Ks_Max))]

							Optimization = BlackBoxOptim.bboptimize(P ->OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; N=P[1], Ψvg=10.0^P[2], Ks=10.0^P[3])[1]; SearchRange = SearchRange, NumDimensions=opt.N_ParamOpt, TraceMode=:silent)
		
							hydro.N[iSoil] 		= BlackBoxOptim.best_candidate(Optimization)[1]
							hydro.Ψvg[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
							hydro.Ks[iSoil] 	= 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[3]

						else
							error( " SoilWater-Toolbox error: option.hydro not found $(option.hydro)")
						end #Option


						Of, Of_θΨ, Of_Kunsat = ofHydro.OF_WRC_KUNSAT(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro) 

						hydro.Nse[iSoil] = 1.0 - Of
						hydro.Nse_θψ[iSoil] = 1.0 - Of_θΨ
						hydro.Nse_Kψ[iSoil] = 1.0 - Of_Kunsat

				end  # for iSoil=1:N_SoilSelect
				
				return hydro
			end  # function: HYDROPARAM_OPT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		vg FUNCTION : HYDRO_PROCESS
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			function OBJECTIVE_FUNCTION(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; N=hydro.N[iSoil], Ψvg=hydro.Ψvg[iSoil], θr=hydro.θr[iSoil], θs=hydro.θs[iSoil], Ks=hydro.Ks[iSoil])

                hydro.θs[iSoil]  = θs
                hydro.θr[iSoil]  = θr
                hydro.Ks[iSoil]  = Ks
                hydro.N[iSoil]   = N
                hydro.Ψvg[iSoil] = Ψvg

				Of, Of_θΨ, Of_Kunsat = ofHydro.OF_WRC_KUNSAT(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro) 
	
				return Of, Of_θΨ, Of_Kunsat
			end  # function: HYDRO_PROCESS
			
		end  # module vg
		# ............................................................
		
end  # module: hydro