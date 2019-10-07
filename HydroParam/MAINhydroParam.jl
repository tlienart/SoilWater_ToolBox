# =============================================================
#		MODULE: name
# =============================================================
module mainHydroParam
	using ..option, ..psdThetar, ..wrc, ..kunsat

	using BlackBoxOptim
	export MAIN_HYDROPARAM

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
	#		FUNCTION : MAIN_HYDROPARAM()
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function MAIN_HYDROPARAM(N_SoilSelect, ∑Psd, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ)

			# WHAT TO OPTIMIZE (just initializing)
				Opt_θs 		= true
				Opt_θr 		= true
				Opt_Ks		= true

			# INITIALIZING
				θs 		= Array{Float64}(undef, (N_SoilSelect))
				θr 		= Array{Float64}(undef, (N_SoilSelect))
				Ks		= Array{Float64}(undef, (N_SoilSelect))
				θsMat 	= Array{Float64}(undef, (N_SoilSelect))

				if option.HydroModel == "Kosugi"
					σ 		= Array{Float64}(undef, (N_SoilSelect))
					Ψm 		= Array{Float64}(undef, (N_SoilSelect))
					σMac 	= Array{Float64}(undef, (N_SoilSelect))
					ΨmMac 	= Array{Float64}(undef, (N_SoilSelect))
					
					hydro = KOSUGI(θs, θr, Ks, σ, Ψm, θsMat, σMac, ΨmMac)

					# Number of parameters to be optimized
						if option.UnimodalBimodal == "Unimodal"
							# Keep constant
							σMac 	= 1.0
							ΨmMac 	= 1000.0
							
							# Parameters to be keept constant
							N_ParamOpt = 5

						elseif option.UnimodalBimodal == "Bimodal"
							# Parameters to be keept constant
							N_ParamOpt = 8

						end  # if: option.UnimodalBimodal == "Unimodal"
					
				elseif option.HydroModel == "Vangenuchten"
					N		= Array{Float64}(undef, (N_SoilSelect))
					Ψvg		= Array{Float64}(undef, (N_SoilSelect))

					hydro = VANGENUCHTEN(θs, θr, Ks, N, Ψvg)

					# Parameters to be optimized
						N_ParamOpt = 5

						Opt_θs 		= true
						Opt_θr 		= true
						Opt_N 		= true
						Opt_Ψvg		= true
						Opt_Ks		= true
				end # option.HydroModel

				θ_Max 		= Array{Float64}(undef, (N_SoilSelect))
				θ_Min 		= Array{Float64}(undef, (N_SoilSelect))
				K_KΨ_Max 	= Array{Float64}(undef, (N_SoilSelect))


			# INITIALIZING
				for iSoil=1:N_SoilSelect
					# LIMITS
						θ_Max[iSoil] 	= maximum(θ_θΨ[iSoil, 1:N_θΨ[iSoil]])  # Greatest measure θ
						θ_Min[iSoil] 	= minimum(θ_θΨ[iSoil, 1:N_θΨ[iSoil]])  # Smallest measure θ
						K_KΨ_Max[iSoil] = maximum(K_KΨ[iSoil, 1:N_KΨ[iSoil]]) # Greatest measure of Kunsat

					# Derive θr from data
						if (option.hydro.θrOpt == "Psd") && (option.Psd)
							hydro.θr[iSoil] = min( psdThetar.PSD_2_θr(iSoil, ∑Psd), θ_Min[iSoil] )
							N_ParamOpt -= 1
							Opt_θr = false
						elseif option.hydro.θrOpt == "Cst"
							# Keep θr = Cst
							θr = param.hydro.θr
							N_ParamOpt -= 1
							Opt_θr = false
						end # option.hydro.θrOpt == "Psd"


					# Derive θs from data
						if option.hydro.θsOpt == "Fixed"
							hydro.θs[iSoil] = θ_Max[iSoil]
							N_ParamOpt -= 1
							Opt_θs = false
						elseif option.hydro.θsOpt == "Φ"
							hydro.θs[iSoil] = param.hydro.Coeff_Φ_2_θs * θ_Max[iSoil]
							Opt_θs = false
						end # option.hydro.θsOpt


					# Derive Ks from data
						if option.hydro.KsOpt == "Fixed" && option.KunsatΨ
							hydro.Ks[iSoil] = K_KΨ_Max[iSoil]
							N_ParamOpt -= 1
							Opt_Ks = false
						elseif !(option.KunsatΨ) # Not optimized
							N_ParamOpt -= 1
							Opt_Ks = false
						end # option.KunsatΨ

					# OPTIMIZING	
						if option.HydroModel == "Kosugi"
							kg.OPTIMISATION_HYDRO(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, Opt_Ks, Opt_θs, Opt_θr, hydro)
						# elseif option.HydroModel == "Vangenuchten"
						# 	vg.OPTIMISATION_HYDRO(N_SoilSelect, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro)
						end

				end  # for iSoil=1:N_SoilSelect


		end  # function: MAIN_HYDROPARAM()


		# =============================================================
		#		MODULE: kg
		# =============================================================
		module kg
			using ...ofHydro
			export OPTIMISATION_HYDRO

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : INITIALIZING_HYDRO
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OPTIMISATION_HYDRO(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, Opt_Ks, Opt_θs, Opt_θr, hydro)


					# Feasible range
					θr_Max[iSoil] = max( min(θ_Min-0.005, param.hydro.θr_Max), 0.0 ) # Maximum value of θr

					θsMat_Min[iSoil] = θ_Max * param.hydro.Coeff_θs_2_θsMat
					θsMat_Max[iSoil] = θ_Max


					SearchRange =[(param.hydro.σ_Min, param.hydro.σ_Max), (param.hydro.Ψm_Min, param.hydro.Ψm_Max), (param.hydro.θr_Min, θr_Max), (param.hydro.θs_Min, param.hydro.θs_Max), (param.hydro.Ks_Min, param.hydro.Ks_Max), (θsMat_Min[iSoil], θsMat_Max[iSoil]), (param.hydro.σMac_Min, param.hydro.σMac_Max), (param.hydro.ΨmMac_Min, param.hydro.ΨmMac_Max)]

					Of = ofHydro.OF_WRC_KUNSAT(iSoil, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro; σ=hydro.σ, Ψm=hydro.Ψm, θr=hydro.θr, θs=hydro.θs, Ks=hydro.Ks, θsMat=hydro.θsMat, σMac=hydro.σMac, ΨmMac=hydro.ΨmMac)
					
					# println("$iSoil, $Of")
		

				return hydro
			end  # function: OPTIMISATION_HYDRO

			
		end  # module kg
		# ............................................................


			# OPTIMIZATION
				# if option.HydroModel == "Kosugi"


				

				# 	Optimization = BlackBoxOptim.bboptimize(P -> ofHydro.OF_WRC_KUNSAT() ; SearchRange =[ (param.θr_Min, θr_Max), (param.σMat_Min, param.σMat_Max), (log10.(param.ΨkgMat_Min), log10.(param.ΨkgMat_Max)) , (θsMat_Min, θsMac[iSoil]), (param.σMac_Min, param.σMac_Max), (log10.(K_Kθ_Min), log10.(param.Ks_Mac_Max))], NumDimensions=N_ParamOpt, TraceMode=:silent)

				# # Values of the optimal bimodal hydraulic params===========================
				# θr[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
				# σMat[iSoil] =  BlackBoxOptim.best_candidate(Optimization)[2]
				# ΨkgMat[iSoil] =  10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[3])
				# θsMat[iSoil] = BlackBoxOptim.best_candidate(Optimization)[4]
				# σMac[iSoil]  = BlackBoxOptim.best_candidate(Optimization)[5]
				# ΨkgMac[iSoil] = param.ΨkgMac

				# if Option_Data_Kθ
			   	# 	KsMac[iSoil] = 10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[6])
				# else
			   	# 	KsMac[iSoil] = -999.0
				# end

				# 	Of_Hydro, hydro = kg.HYDRO_OPTIMISATION()
					

				# elseif option.HydroModel == "Vangenuchten"

				# 	Of_Hydro, hydro = vg.HYDRO_OPTIMISATION()

				# end # option.HydroModel



	# =============================================================
	#		MODULE: vg
	# =============================================================
	module vg
		using ...ofHydro
		export OPTIMISATION_HYDRO

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : INITIALIZING_HYDRO
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function  OPTIMISATION_HYDRO(N_SoilSelect, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro)
			
			
			return hydro
		end  # function: OPTIMISATION_HYDRO
		
	end  # module vg
		# ............................................................
	
end  # module: mainHydroParam
# ............................................................
	
#

# 				θsMat_Min = θsMac[iSoil] * 0.7

# 				Optimization = BlackBoxOptim.bboptimize(P -> of.WRC_KUNSAT(Option_Data_Kθ, iSoil, Ψ_θΨ, θ_θΨ, N_θΨ, K_Kθ, Ψ_Kθ, N_Kθ, θsMac[iSoil], P[1], P[2], 10.0^P[3], P[4], P[5], param.ΨkgMac, 10.0^P[6]) ; SearchRange =[ (param.θr_Min, θr_Max), (param.σMat_Min, param.σMat_Max), (log10.(param.ΨkgMat_Min), log10.(param.ΨkgMat_Max)) , (θsMat_Min, θsMac[iSoil]), (param.σMac_Min, param.σMac_Max), (log10.(K_Kθ_Min), log10.(param.Ks_Mac_Max))], NumDimensions=Dimensions, TraceMode=:silent)

# 				# Values of the optimal bimodal hydraulic params===========================
# 				θr[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
# 				σMat[iSoil] =  BlackBoxOptim.best_candidate(Optimization)[2]
# 				ΨkgMat[iSoil] =  10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[3])
# 				θsMat[iSoil] = BlackBoxOptim.best_candidate(Optimization)[4]
# 				σMac[iSoil]  = BlackBoxOptim.best_candidate(Optimization)[5]
# 				ΨkgMac[iSoil] = param.ΨkgMac

# 				if Option_Data_Kθ
# 			   		KsMac[iSoil] = 10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[6])
# 				else
# 			   		KsMac[iSoil] = -999.0
# 				end


# 			elseif option.BimodalKg == false 
# 		 	# DERIVING KOSUGI UNIMODAL HYDRAULIC PARAMETERS <><><><><><><><><><><><><><><><><><><><><>
# 				if Option_Data_Kθ
# 			   		Dimensions=3 # number of optimized parameters
# 				else
# 			   		Dimensions=4 # number of optimized parameters (one more because of Ks is optimized)
# 				end

# 				θsMat_Min = θsMac[iSoil] * 0.7

# 				σMac_1 = 1.0

# 				Optimization = BlackBoxOptim.bboptimize(P -> of.WRC_KUNSAT(Option_Data_Kθ, iSoil, Ψ_θΨ, θ_θΨ, N_θΨ, K_Kθ, Ψ_Kθ, N_Kθ, θsMac[iSoil], P[1], P[2], 10.0^P[3], θsMac[iSoil], σMac_1, param.ΨkgMac, 10.0^P[4]) ; SearchRange =[ (param.θr_Min, θr_Max), (param.σMat_Min, param.σMat_Max), (log10.(param.ΨkgMat_Min), log10.(param.ΨkgMat_Max)), (log10.(K_Kθ_Min), log10.(param.Ks_Mac_Max))], NumDimensions=Dimensions, TraceMode=:silent)

# 				# Values of the optimal unimodal hydraulic params===========================
# 				θr[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
# 				σMat[iSoil] =  BlackBoxOptim.best_candidate(Optimization)[2]
# 				ΨkgMat[iSoil] =  10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[3])
# 				θsMat[iSoil] = θsMac[iSoil] 
# 				σMac[iSoil]  = σMac_1
# 				ΨkgMac[iSoil] = ΨkgMat[iSoil] 

# 				if Option_Data_Kθ
# 			   		KsMac[iSoil] = 10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[4])
# 				else
# 			   		KsMac[iSoil] = -999.0
# 				end

# 			end # option.BimodalKg


# 			# STATISTICS <><><><><><><><><><><><><><><><><><><><><>
# 				# How good the fit of θ(h)
# 			θ_Obs = zeros(Float64, N_θΨ[iSoil])
# 			θ_Sim_Uni = zeros(Float64, N_θΨ[iSoil])
# 			θ_Sim_Bim = zeros(Float64, N_θΨ[iSoil])
# 			H_Obs = zeros(Float64, N_θΨ[iSoil])

# 			@simd for iH in 1:N_θΨ[iSoil]
# 				H_Obs[iH] = Ψ_θΨ[iSoil,iH]
# 				θ_Obs[iH] = θ_θΨ[iSoil,iH]#

# 				θ_Sim_Uni[iH] = wrc.kg.Ψ_2_θdual(H_Obs[iH], θsMat[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMat[iSoil], σMat[iSoil]) # Unimodal Kosugi WR
# 				θ_Sim_Bim[iH] = wrc.kg.Ψ_2_θdual(H_Obs[iH], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil]) # Bimodal Kosugi WR
# 			end

# 			Nse_θΨ_Uni[iSoil]= 1.0 - stat.NASH_SUTCLIFFE_ERRORmin(θ_Obs[1:N_θΨ[iSoil]], θ_Sim_Uni[1:N_θΨ[iSoil]])
# 			Nse_θΨ_Bim[iSoil]= 1.0 - stat.NASH_SUTCLIFFE_ERRORmin(θ_Obs[1:N_θΨ[iSoil]], θ_Sim_Bim[1:N_θΨ[iSoil]])

# 			# How good the fit of K(θ)
# 			Kunsat_Obs = zeros(Float64,(N_Kθ[iSoil]))
# 			Kunsat_Sim_Uni = zeros(Float64, (N_Kθ[iSoil]))
# 			Kunsat_Sim_Bim = zeros(Float64, (N_Kθ[iSoil]))
# 			@simd for iH in 1:N_Kθ[iSoil]
# 				# Observed
# 				H_Obs = Ψ_Kθ[iSoil,iH]
# 				Kunsat_Obs[iH] = K_Kθ[iSoil,iH]

# 				# Unimodal
# 				Se = wrc.kg.Ψ_2_Se(H_Obs, ΨkgMat[iSoil], σMat[iSoil])
# 				Kunsat_Sim_Uni[iH] = kunsat.kg.Se_2_KUNSAT(Se, θsMat[iSoil], θr[iSoil], σMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMat[iSoil])

# 				# Bimmodal
# 				θ_Sim_Bim = wrc.kg.Ψ_2_θdual(H_Obs, θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil])
# 				Se_Bim = wrc.se.θ_2_Se(θ_Sim_Bim, θsMac[iSoil], θr[iSoil])
# 				Kunsat_Sim_Bim[iH] = kunsat.kg.Se_2_KUNSAT(Se_Bim, θsMac[iSoil], θr[iSoil], σMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMac[iSoil])
# 			end

# 			if Option_Data_Kθ
# 			   Nse_Kθ_Uni[iSoil] = 1.0 - stat.NASH_SUTCLIFFE_ERRORmin(log.(1.0 .+ Kunsat_Obs[1:N_Kθ[iSoil]]), log.(1.0 .+ Kunsat_Sim_Uni[1:N_Kθ[iSoil]]))
# 			   Nse_Kθ_Bim[iSoil] = 1.0 - stat.NASH_SUTCLIFFE_ERRORmin(log.(1.0 .+ Kunsat_Obs[1:N_Kθ[iSoil]]), log.(1.0 .+ Kunsat_Sim_Bim[1:N_Kθ[iSoil]]))
# 			else
# 			   Nse_Kθ_Uni[iSoil] = 0.0
# 			   Nse_Kθ_Bim[iSoil] = 0.0
# 			end
# 		end # looping over soils

# 		return θsMat, θr, σMat, ΨkgMat, KsMat, θsMac, σMac, ΨkgMac, KsMac, Nse_θΨ_Uni, Nse_θΨ_Bim, Nse_Kθ_Uni, Nse_Kθ_Bim 

# 	end  # function OPTIMISE_CHARAC_UNSAT

# end # module optimiseCharacUnsat