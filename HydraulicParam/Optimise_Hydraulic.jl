module optimiseCharacUnsat
	include("Path.jl")
	include("Option.jl")
	include("Cst.jl")
	include("Param.jl")
	include("Reading.jl")
	include("WRC.jl")
	include("Kunsat.jl")
	include("Stats.jl")
	include("ObjectiveFunction.jl")
	include("Table.jl")

	using BlackBoxOptim
	export OPTIMISE_CHARAC_UNSAT

	 # =========================================
	 #        OPTIMISE_CHARAC_UNSAT
	 # =========================================
   function OPTIMISE_CHARAC_UNSAT(Ψ_θΨ, θ_θΨ, N_θΨ, Nsample, Φ; N_Kθ=ones(Int8, Nsample), K_Kθ=zeros(Float64,Nsample,1), Ψ_Kθ=zeros(Float64,Nsample,1), θr_Psd=zeros(Float64, Nsample), Option_Data_Kθ = true )
	"""
		Φorθs = "θs" OR "Φ"
			IF Φorθs = "θs" than Φ = θs
			IF Φorθs = "Φ" than Φ is computed
	"""

		# Putting in memory
		σMat = zeros(Float64, Nsample)
		σMac = zeros(Float64, Nsample)
		θsMat = zeros(Float64, Nsample)
		θsMac = zeros(Float64, Nsample)
		θr = zeros(Float64, Nsample)
		Nse_θΨ_Uni = zeros(Float64, Nsample)
		Nse_θΨ_Bim = zeros(Float64, Nsample)
		Nse_Kθ_Uni = zeros(Float64, Nsample)
		Nse_Kθ_Bim = zeros(Float64, Nsample)
		KsMac = zeros(Float64, Nsample)
		KsMat = zeros(Float64, Nsample)
		ΨkgMat = zeros(Float64, Nsample)
		ΨkgMac = zeros(Float64, Nsample)

		@simd for iSoil in 1:Nsample

			# ================= LIMITS =================
				θ_Max = maximum(θ_θΨ[iSoil,1:N_θΨ[iSoil]])  # Greatest measured θ before Φ
				θ_Min = minimum(θ_θΨ[iSoil,1:N_θΨ[iSoil]])  # Smallest measured θ

				if option.Φorθs == "Φ"
					θsMac[iSoil] = max( param.Coeff_Φ_2_θs * Φ[iSoil], θ_Max + 0.005 ) # Deriving θsMac
				else
					θsMac[iSoil] = Φ[iSoil]
				end

				θr_Max = max( min(θ_Min-0.005, param.θr_Max), 0. ) # Maximum value of θr
				K_Kθ_Min = minimum(K_Kθ[iSoil,1:N_Kθ[iSoil]]) # Greatest measured k(h) before KsMac
			# ======================================

		 # DERIVING KOSUGI UNIMODAL HYDRAULIC PARAMETERS <><><><><><><><><><><><><><><><><><><><><>
			if Option_Data_Kθ
			   Dimensions=5
			else
			   Dimensions=4
			end

		 	if option.BimodalKg == true	
			 # DERIVING KOSUGI BIMODAL HYDRAULIC PARAMETERS <><><><><><><><><><><><><><><><><><><><><>
				if Option_Data_Kθ
			  	 	Dimensions=5 # number of optimized parameters
				else
			   		Dimensions=6 # number of optimized parameters (one more because of Ks is optimized)
				end

				θsMat_Min = θsMac[iSoil] * 0.7

				Optimization = BlackBoxOptim.bboptimize(P -> of.WRC_KUNSAT(Option_Data_Kθ, iSoil, Ψ_θΨ, θ_θΨ, N_θΨ, K_Kθ, Ψ_Kθ, N_Kθ, θsMac[iSoil], P[1], P[2], 10.0^P[3], P[4], P[5], param.ΨkgMac, 10.0^P[6]) ; SearchRange =[ (param.θr_Min, θr_Max), (param.σMat_Min, param.σMat_Max), (log10.(param.ΨkgMat_Min), log10.(param.ΨkgMat_Max)) , (θsMat_Min, θsMac[iSoil]), (param.σMac_Min, param.σMac_Max), (log10.(K_Kθ_Min), log10.(param.Ks_Mac_Max))], NumDimensions=Dimensions, TraceMode=:silent)

				# Values of the optimal bimodal hydraulic params===========================
				θr[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
				σMat[iSoil] =  BlackBoxOptim.best_candidate(Optimization)[2]
				ΨkgMat[iSoil] =  10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[3])
				θsMat[iSoil] = BlackBoxOptim.best_candidate(Optimization)[4]
				σMac[iSoil]  = BlackBoxOptim.best_candidate(Optimization)[5]
				ΨkgMac[iSoil] = param.ΨkgMac

				if Option_Data_Kθ
			   		KsMac[iSoil] = 10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[6])
				else
			   		KsMac[iSoil] = -999.0
				end


			elseif option.BimodalKg == false 
		 	# DERIVING KOSUGI UNIMODAL HYDRAULIC PARAMETERS <><><><><><><><><><><><><><><><><><><><><>
				if Option_Data_Kθ
			   		Dimensions=3 # number of optimized parameters
				else
			   		Dimensions=4 # number of optimized parameters (one more because of Ks is optimized)
				end

				θsMat_Min = θsMac[iSoil] * 0.7

				σMac_1 = 1.0

				Optimization = BlackBoxOptim.bboptimize(P -> of.WRC_KUNSAT(Option_Data_Kθ, iSoil, Ψ_θΨ, θ_θΨ, N_θΨ, K_Kθ, Ψ_Kθ, N_Kθ, θsMac[iSoil], P[1], P[2], 10.0^P[3], θsMac[iSoil], σMac_1, param.ΨkgMac, 10.0^P[4]) ; SearchRange =[ (param.θr_Min, θr_Max), (param.σMat_Min, param.σMat_Max), (log10.(param.ΨkgMat_Min), log10.(param.ΨkgMat_Max)), (log10.(K_Kθ_Min), log10.(param.Ks_Mac_Max))], NumDimensions=Dimensions, TraceMode=:silent)

				# Values of the optimal unimodal hydraulic params===========================
				θr[iSoil] = BlackBoxOptim.best_candidate(Optimization)[1]
				σMat[iSoil] =  BlackBoxOptim.best_candidate(Optimization)[2]
				ΨkgMat[iSoil] =  10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[3])
				θsMat[iSoil] = θsMac[iSoil] 
				σMac[iSoil]  = σMac_1
				ΨkgMac[iSoil] = ΨkgMat[iSoil] 

				if Option_Data_Kθ
			   		KsMac[iSoil] = 10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[4])
				else
			   		KsMac[iSoil] = -999.0
				end

			end # option.BimodalKg


			# STATISTICS <><><><><><><><><><><><><><><><><><><><><>
				# How good the fit of θ(h)
			θ_Obs = zeros(Float64, N_θΨ[iSoil])
			θ_Sim_Uni = zeros(Float64, N_θΨ[iSoil])
			θ_Sim_Bim = zeros(Float64, N_θΨ[iSoil])
			H_Obs = zeros(Float64, N_θΨ[iSoil])

			@simd for iH in 1:N_θΨ[iSoil]
				H_Obs[iH] = Ψ_θΨ[iSoil,iH]
				θ_Obs[iH] = θ_θΨ[iSoil,iH]#

				θ_Sim_Uni[iH] = wrc.kg.Ψ_2_θdual(H_Obs[iH], θsMat[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMat[iSoil], σMat[iSoil]) # Unimodal Kosugi WR
				θ_Sim_Bim[iH] = wrc.kg.Ψ_2_θdual(H_Obs[iH], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil]) # Bimodal Kosugi WR
			end

			Nse_θΨ_Uni[iSoil]= 1.0 - stats.NASH_SUTCLIFFE_ERRORmin(θ_Obs[1:N_θΨ[iSoil]], θ_Sim_Uni[1:N_θΨ[iSoil]])
			Nse_θΨ_Bim[iSoil]= 1.0 - stats.NASH_SUTCLIFFE_ERRORmin(θ_Obs[1:N_θΨ[iSoil]], θ_Sim_Bim[1:N_θΨ[iSoil]])

			# How good the fit of K(θ)
			Kunsat_Obs = zeros(Float64,(N_Kθ[iSoil]))
			Kunsat_Sim_Uni = zeros(Float64, (N_Kθ[iSoil]))
			Kunsat_Sim_Bim = zeros(Float64, (N_Kθ[iSoil]))
			@simd for iH in 1:N_Kθ[iSoil]
				# Observed
				H_Obs = Ψ_Kθ[iSoil,iH]
				Kunsat_Obs[iH] = K_Kθ[iSoil,iH]

				# Unimodal
				Se = wrc.kg.Ψ_2_Se(H_Obs, ΨkgMat[iSoil], σMat[iSoil])
				Kunsat_Sim_Uni[iH] = kunsat.kg.Se_2_KUNSAT(Se, θsMat[iSoil], θr[iSoil], σMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMat[iSoil])

				# Bimmodal
				θ_Sim_Bim = wrc.kg.Ψ_2_θdual(H_Obs, θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil])
				Se_Bim = wrc.se.θ_2_Se(θ_Sim_Bim, θsMac[iSoil], θr[iSoil])
				Kunsat_Sim_Bim[iH] = kunsat.kg.Se_2_KUNSAT(Se_Bim, θsMac[iSoil], θr[iSoil], σMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMac[iSoil])
			end

			if Option_Data_Kθ
			   Nse_Kθ_Uni[iSoil] = 1.0 - stats.NASH_SUTCLIFFE_ERRORmin(log.(1.0 .+ Kunsat_Obs[1:N_Kθ[iSoil]]), log.(1.0 .+ Kunsat_Sim_Uni[1:N_Kθ[iSoil]]))
			   Nse_Kθ_Bim[iSoil] = 1.0 - stats.NASH_SUTCLIFFE_ERRORmin(log.(1.0 .+ Kunsat_Obs[1:N_Kθ[iSoil]]), log.(1.0 .+ Kunsat_Sim_Bim[1:N_Kθ[iSoil]]))
			else
			   Nse_Kθ_Uni[iSoil] = 0.0
			   Nse_Kθ_Bim[iSoil] = 0.0
			end
		end # looping over soils

		return θsMat, θr, σMat, ΨkgMat, KsMat, θsMac, σMac, ΨkgMac, KsMac, Nse_θΨ_Uni, Nse_θΨ_Bim, Nse_Kθ_Uni, Nse_Kθ_Bim 

	end  # function OPTIMISE_CHARAC_UNSAT

end # module optimiseCharacUnsat