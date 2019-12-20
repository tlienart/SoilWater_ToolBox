module quasiExact # quasi-exact objective function
	import ..option, ..sorptivity, ..wrc, ..kunsat, ..option, ..param
	import BlackBoxOptim, Optim
 	export INFILTRATION3D_2_1D, OF_INFILTRATION_2_HYDRO

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_2_TIMEη
	#		= COMPUTE NORMALISED TIME: Tη =
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_2_TIMEη(Sorptivity, T, ΔK)
			return Time_η = T * 2.0 * (ΔK / Sorptivity) ^ 2.0
		end # function: TIME_2_TIMEη


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION3D_2_1D
	# 		= TRANSFORMS INFILTRATION_3D TO INFILTRATION_1D =
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION3D_2_1D(hydroInfilt, ∑Infilt, infiltParam, iSoil, Sorptivity, T; ϵ=eps())
			Δθ = hydroInfilt.θs[iSoil] - hydroInfilt.θr[iSoil]
			return ∑Infilt_1D = max(∑Infilt - (T * infiltParam.γ[iSoil] * Sorptivity ^ 2.0) / (infiltParam.RingRadius[iSoil] * Δθ), ϵ)
		end # function : INFILTRATION3D_2_1D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  INFILTTRATIONη_2_3D
	# 		TRANSFORMS NORMALIZED INFILTRATION TO INFILTRATION-3D
	#		Function compute infiltration-1d from normalized infiltration
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTTRATIONη_2_3D(Infilt_η, infiltParam, iSoil, K_θini, Sorptivity, T, Time_η, ΔK, Δθ)
			
			INFILTTRATIONη_2_1D(Infilt_η, K_θini, Sorptivity, T, ΔK) = K_θini * T + Infilt_η * (Sorptivity ^ 2.0) / (2.0 * ΔK)

			ΔI_η = Time_η * infiltParam.γ[iSoil]

			return ∑Infilt = INFILTTRATIONη_2_1D(Infilt_η, K_θini, Sorptivity, T, ΔK) + ΔI_η * (Sorptivity ^ 4.0) / (2.0 * infiltParam.RingRadius[iSoil] * Δθ * (ΔK ^ 2.0))
		end # function :  INFILTTRATIONη_2_3D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION3D_2_η
 	# 		TRANSFORMS INFILTRATION-3D TO NORMALIZED INFILTRATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION3D_2_η(∑Infilt, infiltParam, iSoil, K_θini, Sorptivity, T, ΔK, Δθ; ϵ=eps())
			return Infilt_η = max((2.0 * ΔK / Sorptivity ^ 2.0) * (∑Infilt - K_θini * T - infiltParam.γ[iSoil] * TIME_2_TIMEη(Sorptivity, T, ΔK) * (Sorptivity^4.0) / (infiltParam.RingRadius[iSoil] * Δθ * 2.0* (ΔK^2.0)) ), ϵ)
		end # function INFILTRATION3D_2_η


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDRO_2_INFILTRATION3D
	# 		COMPUTE INFILTRATION_3D FROM OPTIMIZED HYDRAULIC PARAMETERS
	# 		Solving quasiexact solution
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
		function HYDRO_2_INFILTRATION3D(∑Infilt, hydroInfilt, infiltParam, iSoil, N_Infilt, T)

			Infilt_η = Array{Float64}(undef, N_Infilt[iSoil])

			Sorptivity = sorptivity.SORPTIVITY(infiltParam.θ_Ini[iSoil], iSoil, hydroInfilt)

			Se_Ini = wrc.θ_2_Se(infiltParam.θ_Ini[iSoil], iSoil, hydroInfilt)

			K_θini = kunsat.Se_2_KUNSAT(Se_Ini, iSoil, hydroInfilt)

			ΔK = hydroInfilt.Ks[iSoil] - K_θini

			Δθ = hydroInfilt.θs[iSoil] - infiltParam.θ_Ini[iSoil]
		
			# At t=1
				∑Infilt[1] = 0.0
				Infilt_η[1] = 0.0
				Infilt_η_Min = 0.0001 # 0.0001
				Infilt_η_Max = 1.0 * T[2] #Since T[1] = 0 0.05

			# ~~~~~~~~~~~~~~~~~~~~
			function OF_QUASIEXACTη(Infilt_η, infiltParam, iSoil, Time_η)
				Left_Term = Time_η

				Right_Term = (1.0 / (1.0 - infiltParam.β[iSoil])) * (Infilt_η - log((exp(infiltParam.β[iSoil] * Infilt_η) + infiltParam.β[iSoil] - 1.0) / infiltParam.β[iSoil]))
				
				if Right_Term < 0.0
					return OF = 10000.0 * exp(Infilt_η)
				else
					return OF = abs(Left_Term - Right_Term)
				end
			end # function OF_QUASIEXACTη ~~~~~~~~~~~~~~~~~~~~


			for iT in 2:N_Infilt[iSoil] # Looping for every time step
				Time_η = TIME_2_TIMEη(Sorptivity, T[iSoil,iT], ΔK)

				# Solving for Infilt_η
					Optimization = Optim.optimize(Infilt_η -> OF_QUASIEXACTη(Infilt_η, infiltParam, iSoil, Time_η), Infilt_η_Min, Infilt_η_Max, Optim.GoldenSection())
					Infilt_η[iT] = Optim.minimizer(Optimization)[1]

				# Deriving the new bounds such that infiltration increases with time & the slope decreases with time
					Infilt_η_Min = Infilt_η[iT] + 0.0001

				# Maximum infiltration rate for T+1: (Infilt[T2] - Infilt[T1]) / (T2 - T1) which is 1 seconds
					if iT <= N_Infilt[iSoil] - 1
						Infilt_η_Max = Infilt_η[iT] + (T[iSoil,iT+1]- T[iSoil,iT]) * (Infilt_η[iT] - Infilt_η[iT-1]) / (T[iSoil,iT] - T[iSoil,iT-1])
					else
						Infilt_η_Max = Infilt_η[iT] + (Infilt_η[iT] - Infilt_η[iT-1])
						
					end

				# Transforming INFILTTRATIONη to INFILTRATION3D 
					∑Infilt[iSoil,iT] =  INFILTTRATIONη_2_3D(Infilt_η[iT], infiltParam, iSoil, K_θini, Sorptivity, T[iSoil,iT], Time_η, ΔK, Δθ)
			end
			return ∑Infilt
		end # function: HYDRO_2_INFILTRATION3D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_INFILTRATION_2_HYDRO
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_INFILTRATION_2_HYDRO(∑Infilt_Obs, infiltOutput, infiltParam, iSoil, N_Infilt, T, hydroInfilt; σ=hydroInfilt.σ[iSoil], Ψm=hydroInfilt.Ψm[iSoil], θr=hydroInfilt.θr[iSoil], θs=hydroInfilt.θs[iSoil], Ks=hydroInfilt.Ks[iSoil], W=0.9)

			hydroInfilt.θs[iSoil] = θs
			hydroInfilt.θr[iSoil] = θr
			hydroInfilt.Ks[iSoil] = Ks
			hydroInfilt.σ[iSoil] = σ
			hydroInfilt.Ψm[iSoil] = Ψm
			hydroInfilt.θsMat[iSoil] = θs
			hydroInfilt.ΨmMac[iSoil] = Ψm
			hydroInfilt.σMac[iSoil] = σ

			iT_TransSteady = infiltOutput.iT_TransSteady_Data[iSoil]

			Sorptivity = sorptivity.SORPTIVITY(infiltParam.θ_Ini[iSoil], iSoil, hydroInfilt)

			Se_Ini = wrc.θ_2_Se(infiltParam.θ_Ini[iSoil], iSoil, hydroInfilt)

			K_θini = kunsat.Se_2_KUNSAT(Se_Ini, iSoil, hydroInfilt)

			ΔK = hydroInfilt.Ks[iSoil] - K_θini

			Δθ = hydroInfilt.θs[iSoil] - infiltParam.θ_Ini[iSoil]

			Left_Term = zeros(Float64, N_Infilt[iSoil])
			Right_Term = zeros(Float64, N_Infilt[iSoil])
			Of_Penalty = 0.0 ; Of_Stead = 0.0 ; Of_Trans = 0.0
			for iT in 2:N_Infilt[iSoil]
				Time_η = TIME_2_TIMEη(Sorptivity, T[iSoil,iT], ΔK)

				Infilt_η = INFILTRATION3D_2_η(∑Infilt_Obs[iSoil,iT], infiltParam, iSoil, K_θini, Sorptivity, T[iSoil,iT], ΔK, Δθ)

				Left_Term[iT] = Time_η

				Right_Term[iT] = (1.0 / (1.0 - infiltParam.β[iSoil])) * (Infilt_η - log((exp(infiltParam.β[iSoil] * Infilt_η) + infiltParam.β[iSoil] - 1.0) / infiltParam.β[iSoil]))
				
				if Right_Term[iT] > 0.0
					if iT <= infiltOutput.iT_TransSteady_Data[iSoil]
						Of_Trans += ((Left_Term[iT]) - (Right_Term[iT])) ^ 2.0
					else
						Of_Stead += (log10(Left_Term[iT]) - log10(Right_Term[iT])) ^ 2.0
					end #  iT <= infiltOutput.iT_TransSteady_Data
				else
					Of_Penalty += 1000.0 * exp(Infilt_η)
					Right_Term[iT] = 0.0
				end #  Right_Term[iT] > 0.0
			end #  Right_Term[iT] > 0.0

			return Wof = (W * Of_Trans / Float64(iT_TransSteady-1)) + ((1.0 - W) * Of_Stead / Float64(N_Infilt[iSoil] - iT_TransSteady + 1)) + Of_Penalty
		end # function: OF_INFILTRATION_2_HYDRO


		# #= =============== KOSUGI =============== =#
	# module kg	
	# 	import ..option, ..sorptivity, ..wrc, ..kunsat, ..param
	# 	import ..quasiExact
	# 	export BlackBoxOptim, Optim
	# 	export INFILTRATION3D_2_HYDRO, INFILTRATION3D_2_HYDRO_σMOD

	# 	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# 	#		FUNCTION : INFILTRATION3D_2_HYDRO_σMOD
	# 	# 		OPTIMIZE THE HYDRAULIC PARAMETERS FROM INFILTRATION DATA
	# 	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	# 		# function INFILTRATION3D_2_HYDRO_σMOD(T,  ∑Infilt_Obs, N_Infilt, θs, θ_Ini, Time_TransStead, Sorptivity, infiltParam, hydroInfilt)
	# 		# 	θr = 0.0
	# 		# 	Δθ = θs - θ_Ini
	# 		# 	Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)

	# 		# 	# Objective function which matches observed with simulated infiltration
	# 		# 	function OF_Fast_η(T,  ∑Infilt_Obs, Δθ, N_Infilt, infiltParam, θr, θ_Ini, σ, iT_TransStead)
	# 		# 		Hkg_σ = relationship.σ_2_Hkg(σ)

	# 		# 		Hkg_Sorpt = sorptivity.kg.SORPTIVITY_2_Hkg(Sorptivity, θ_Ini, θs, θr, σ, Ks)

	# 		# 		K_θini = kunsat.kg.Se_2_KUNSAT(Se_Ini, θs, θr, σ, Ks, θs, σ)

	# 		# 		ΔK = Ks - K_θini

	# 		# 		Of_Hkg = abs(log10(Hkg_Sorpt) - log10(Hkg_σ))

	# 		# 		Wof = quasiExact.OF_INFILTRATION_2_HYDRO( ∑Infilt_Obs[1:N_Infilt[iSoil]], infiltOutput, infiltParam, iSoil, K_θini, N_Infilt, Sorptivity, T[1:N_Infilt[iSoil]], ΔK, Δθ) + Of_Hkg / 10.0

	# 		# 		return Wof
	# 		# 	end # function INFILTRATION3D_2_HYDRO_σMOD

	# 		# 	# OPTIMIZATION

	# 		# 		Optimization = Optim.optimize(σ -> OF_Fast_η(T[1:N_Infilt[iSoil]],  ∑Infilt_Obs[1:N_Infilt[iSoil] ], Δθ, N_Infilt , infiltParam, θr, θ_Ini, σ, iT_TransStead), param.σ_Min, param.σ_Max, GoldenSection() )

	# 		# 		# Values of the optimal hydraulic params
	# 		# 		σ = Optim.minimizer(Optimization)[1]
				
			
	# 		# 	Hkg_Sorpt = sorptivity.kg.SORPTIVITY_2_Hkg(Sorptivity, θ_Ini, θs, θr, σ, Ks)
			
	# 		# 	return σ, Hkg_Sorpt
	# 		# end # function INFILTRATION3D_2_HYDRO_σMOD





	# end # module kg


	#= =============== VAN GENUCHTEN =============== =#
	# module vg
	# include("Cst.jl")
	# include("Param.jl")
	# include("Wrc.jl")
	# include("Kunsat.jl")
	# include("Sorptivity.jl")
	# include("Stat.jl")
	# include("Tools.jl")
	# include("HydroRelationship.jl")
	# using ..quasiExact
	# using BlackBoxOptim, Optim
	# 	export INFILTRATION3D_2_HYDRO

	# 	# OPTIMIZE THE HYDRAULIC PARAMETERS FROM INFILTRATION DATA============================
	# 	function INFILTRATION3D_2_HYDRO(T,  ∑Infilt_Obs, N_Infilt, θs, θ_Ini, infiltParam, Km, iT_TransStead, Time_TransStead, Option_Opt_N, N=1.)
	# 		θr = 0.
	# 		Δθ = θs - θ_Ini
			
	# 		# Objective function which matches observed with simulated infiltration
	# 		function OF_Fast_η(T,  ∑Infilt_Obs, Δθ, N_Infilt, infiltParam, θr, θ_Ini, Hvg, Ks, N, Km)
			
	# 		 Sorptivity = sorptivity.vg.SORPTIVITY(θ_Ini, θs, θs, θr, Hvg, N, Ks, Km)

	# 			Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)
			
	# 			K_θini = kunsat.vg.Se_2_KUNSAT(Se_Ini, N, Ks, Km)
			
	# 			ΔK = Ks - K_θini

	# 			OF_Cumul = quasiExact.OF_INFILTRATION_2_HYDRO(N_Infilt[iSoil], iT_TransStead, T[1:N_Infilt[iSoil]],  ∑Infilt_Obs[1:N_Infilt[iSoil]],Sorptivity, ΔK, K_θini, Δθ, infiltParam)
	# 			return OF_Cumul
	# 		end

	# 		# OPTIMIZATION
	# 		# N_Min and N_Max depends on Km
	# 		if Km == 1
	# 			N_Min = param.N_Km1_Min
	# 			N_Max = param.N_Km1_Max
	# 		elseif Km == 2
	# 			N_Min = param.N_Km2_Min
	# 			N_Max = param.N_Km2_Max
	# 		end

	# 		if Option_Opt_N # If Ks is not known
	# 			Optimization = BlackBoxOptim.bboptimize(Param -> OF_Fast_η(T[1:N_Infilt[iSoil]],  ∑Infilt_Obs[1:N_Infilt[iSoil]], Δθ, N_Infilt[iSoil], infiltParam, θr, θ_Ini, 10.0 ^Param[1], 10.0 ^Param[2], Param[3], Km) ; SearchRange =[ (log10(param.Hvg_Min), log10(param.Hvg_Max)), (log10(param.Ks_Min), log10(param.Ks_Max)), (N_Min, N_Max)], NumDimensions=3, TraceMode=:silent)
	# 			# Values of the optimal hydraulic params
	# 			Hvg = 10.0 ^(BlackBoxOptim.best_candidate(Optimization)[1])
	# 			Ks = 10.0 ^(BlackBoxOptim.best_candidate(Optimization)[2])
	# 			N = BlackBoxOptim.best_candidate(Optimization)[3]
	# 		else # if N is known
	# 			Optimization = BlackBoxOptim.bboptimize(Param -> OF_Fast_η(T[1:N_Infilt[iSoil]],  ∑Infilt_Obs[1:N_Infilt[iSoil]], Δθ, N_Infilt[iSoil], infiltParam, θr, θ_Ini, 10.0 ^Param[1], 10.0 ^Param[2], N, Km) ; SearchRange = [ (log10(param.Hvg_Min), log10(param.Hvg_Max)), (param.Ks_Log_Min, param.Ks_Log_Max)], NumDimensions=2, TraceMode=:silent)
	# 			# Values of the optimal hydraulic params
	# 			Hvg = 10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[1])
	# 			Ks = 10.0 ^(BlackBoxOptim.best_candidate(Optimization)[2])  
	# 		end

	# 	 Sorptivity = sorptivity.vg.SORPTIVITY(θ_Ini, θs, θs, θr, Hvg, N, Ks, Km)

	# 		Se_Ini= wrc.se.θ_2_Se(θ_Ini, θs, θr)
 
	# 		Kr_θini= kunsat.vg.Se_2_KUNSAT(Se_Ini, N, 1., Km)
			
	# 		return Ks, Kr_θini,Sorptivity, N, Hvg
	# 	end
	# end # module vg

end # module quasiExact