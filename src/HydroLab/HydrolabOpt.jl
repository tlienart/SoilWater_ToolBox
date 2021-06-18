# =============================================================
#		module: hypixOpt
# =============================================================
module hydrolabOpt
	import ..ofHydrolab, ..tool, ..optimize, ..hydroRelation, ..psdThetar
	using BlackBoxOptim, Statistics
	export HYDROLABOPT_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIXOPT_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYDROLABOPT_START(;N_iZ, ∑Psd, θ_θΨobs, Ψ_θΨobs, N_θΨobs, K_KΨobs=[0], Ψ_KΨobs=[0], N_KΨobs=1, hydro, hydroOther, option, optionₘ, optim, param, θϵ=0.005)
		for iZ = 1:N_iZ
			# CORRECTION OF THE FEASIBLE RANGE ~~~
				θobs_Min = minimum(θ_θΨobs[iZ, 1:N_θΨobs[iZ]])  	# Smallest measure θ

				θobs_Max = maximum(θ_θΨobs[iZ, 1:N_θΨobs[iZ]])  	# Greatest measure θ

			# CORRECTING Θr ~~~
				if ("θr" ∈ optim.ParamOpt)
					hydro.θr_Max[iZ] = max( min(θobs_Min-θϵ, hydro.θr_Max[iZ]), hydro.θr_Min[iZ] ) # Maximum value of θr

					# Changing the feasible range of θr
					iθr = findfirst(isequal("θr"), optim.ParamOpt)[1]
					optim.ParamOpt_Max[iθr] = hydro.θr_Max[iZ]

				elseif ("θr" ∉ optim.ParamOpt) && (optionₘ.θrOpt==:ParamPsd) && (option.data.Psd) # Derive θr frpm PSD
					hydro.θr[iZ] = min(psdThetar.PSD_2_θr_FUNC(∑Psd, hydro, iZ, param), θobs_Min-θϵ)

				end # if ("θr" ∈ optim.ParamOpt)

			# TEST IF EXIST Ψ=0  ~~~
				if minimum(Ψ_θΨobs[iZ,1:N_θΨobs[iZ]]) < eps(1000.0)
					Flag_Ψ0 = true
				else
					Flag_Ψ0 = false
				end

			# CORRECTING θS  ~~~
				if ("θs" ∈ optim.ParamOpt) && Flag_Ψ0
					hydro.θs_Min[iZ] = θobs_Max * 0.75
					hydro.θs_Max[iZ] = θobs_Max * 1.1
					hydro.Φ[iZ] = θobs_Max / param.hydro.Coeff_Φ_2_θs

					# Changing the feasible range of θs
						iθs = findfirst(isequal("θs"), optim.ParamOpt)[1]
						optim.ParamOpt_Min[iθs] = hydro.θs_Min[iZ]
						optim.ParamOpt_Max[iθs] = hydro.θs_Max[iZ]

				elseif ("θs" ∉ optim.ParamOpt) && Flag_Ψ0 # <>=<>=<>=<>=<>
						hydro.θs[iZ] = θobs_Max
						hydro.Φ[iZ] = hydro.θs[iZ] / param.hydro.Coeff_Φ_2_θs

				elseif  ("θs" ∉ optim.ParamOpt) && (optionₘ.θsOpt==:Φ) # <>=<>=<>=<>=<>
						if hydro.Φ[iZ] * 0.95 > θobs_Max + θϵ
							hydro.θs[iZ] = hydro.Φ[iZ] * 0.95
						elseif hydro.Φ[iZ] * 0.965 > θobs_Max + θϵ
							hydro.θs[iZ] = hydro.Φ[iZ] * 0.965
						else
							hydro.θs[iZ] = max(hydro.Φ[iZ] - θϵ, θobs_Max + θϵ)
						end # hydro.Φ[iZ] * 0.95 > θobs_Max + θϵ

				elseif ("θs" ∈ optim.ParamOpt) && (optionₘ.θsOpt==:FromData) # <>=<>=<>=<>=<>
					hydro.θs_Min[iZ] = θobs_Max - θϵ
					hydro.θs_Max[iZ] = max(θobs_Max + 0.10, hydro.θs_Max[iZ])

					# Changing the feasible range of θs
						iθs = findfirst(isequal("θs"), optim.ParamOpt)[1]
						optim.ParamOpt_Min[iθs] = hydro.θs_Min[iZ]
						optim.ParamOpt_Max[iθs] = hydro.θs_Max[iZ]
		
				end # optionₘ.θsOpt
				
			# CORRECTING Ks  ~~~
				if optionₘ.KsOpt
					# test if exist Ψ=0
					if "Ks" ∈ optim.ParamOpt
						if minimum(Ψ_KΨobs[iZ,1:N_θΨobs[iZ]]) < eps(100.0)
							Flag_K0 = true
						else
							Flag_K0 = false
						end
					end # if "Ks" ∈ optim.ParamOpt

					K_KΨobs_Max = maximum(K_KΨobs[iZ, 1:N_KΨobs[iZ]])

					if !(Flag_K0) && ("Ks" ∈ optim.ParamOpt)
						hydro.Ks_Min[iZ] = K_KΨobs_Max # Greatest measure of Kunsat)

						# Modifying the searchrange
						iKs = findfirst(isequal("Ks"), optim.ParamOpt)[1]
						optim.ParamOpt_Min[iKs] = K_KΨobs_Max

					elseif Flag_K0 && ("Ks" ∈ optim.ParamOpt)
						hydro.Ks_Max[iZ] = K_KΨobs_Max # Greatest measure of Kunsat

						# Modifying the searchrange
							iKs = findfirst(isequal("Ks"), optim.ParamOpt)[1]
							optim.ParamOpt_Max[iKs] = hydro.Ks_Max[iZ]

					elseif ("Ks" ∉ optim.ParamOpt)
						hydro.Ks_Max[iZ]  = K_KΨobs_Max
						hydro.Ks[iZ] = hydro.Ks_Max[iZ]

					end # "Ks" ∈ optim.ParamOpt
				end # if optionₘ.KsOpt
			
			# Updated searchrange
				SearchRange = optimize.SEARCHRANGE(optionₘ, optim)

			# OPTIMIZATION: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				Optimization = BlackBoxOptim.bboptimize(X -> hydrolabOpt.OF_HYDROLAB(X, hydro, iZ, K_KΨobs, N_KΨobs, N_θΨobs, optim, optionₘ, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs); SearchRange=SearchRange, NumDimensions=optim.NparamOpt, TraceMode=:silent)

				X = BlackBoxOptim.best_candidate(Optimization)

				hydro = hydrolabOpt.PARAM_2_hydro(hydro, iZ, optim, optionₘ, X)

				# STATISTICS
					Of, Of_θΨ, Of_Kunsat = ofHydrolab.OF_WRC_KUNSAT(optionₘ, iZ, θ_θΨobs, Ψ_θΨobs, N_θΨobs, K_KΨobs, Ψ_KΨobs, N_KΨobs, hydro, optim) 

					hydroOther.Rmse[iZ], hydroOther.Rmse_KΨ[iZ], hydroOther.Rmse_θΨ[iZ] = ofHydrolab.OF_RMSE(optionₘ, iZ, θ_θΨobs, Ψ_θΨobs, N_θΨobs, K_KΨobs, Ψ_KΨobs, N_KΨobs, hydro, optim) 

					hydroOther.Nse[iZ]    = 1.0 - Of
					hydroOther.Nse_θΨ[iZ] = 1.0 - Of_θΨ

					if optionₘ.KsOpt
						hydroOther.Nse_KΨ[iZ] = 1.0 - Of_Kunsat
					end
		end # for iZ = 1:N_iZ

		# OVERALL STATISTICS OF THE OPTIMIZATION
			Nse_θΨ_Aver = Statistics.mean(hydroOther.Nse_θΨ[1:N_iZ])
			Nse_KΨ_Aver = Statistics.mean(max.(hydroOther.Nse_KΨ[1:N_iZ],0.0))

			Rmse_Aver    = Statistics.mean(hydroOther.Rmse[1:N_iZ])
			Rmse_θΨ_Aver = Statistics.mean(hydroOther.Rmse_θΨ[1:N_iZ])
			Rmse_KΨ_Aver = Statistics.mean(hydroOther.Rmse_KΨ[1:N_iZ])
				
			if optionₘ.KsOpt
				Nse_Aver = (Nse_θΨ_Aver + Nse_KΨ_Aver) / 2.0
			else
				Nse_Aver = Nse_θΨ_Aver
			end

			println("	=== === Optimizing Hydraulic parameters === ")
			println("    		~  Nse_θΨ = $(round(Nse_θΨ_Aver,digits=3)),  Nse_KΨ = $(round(Nse_KΨ_Aver,digits=3)), Nse = $(round(Nse_Aver,digits=3))  ~")
			println("    		~  Rmse_θΨ = $(round(Rmse_θΨ_Aver,digits=4)),  Rlmse_KΨ = $(round(Rmse_KΨ_Aver,digits=4)), Rmse = $(round(Rmse_Aver,digits=4))  ~ \n")
			println( "	=== === ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ === ===")
	return hydro, hydroOther
	end  # function: HYPIXOPT_START


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_HYDROLAB(X, hydro, iZ, K_KΨobs, N_KΨobs, N_θΨobs, optim, optionₘ, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs)
			# New optimized which are put into the matching veg or hydro parameters
				hydro = hydrolabOpt.PARAM_2_hydro(hydro, iZ, optim, optionₘ, X)
		
			# Weighted Objective Function
				Of, Of_θΨ, Of_Kunsat = ofHydrolab.OF_WRC_KUNSAT(optionₘ, iZ, θ_θΨobs, Ψ_θΨobs, N_θΨobs, K_KΨobs, Ψ_KΨobs, N_KΨobs, hydro, optim) 
				
		return Of
		end  # function: OF_HYPIX


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PARAM_2_hydro(hydro, iZ, optim, optionₘ, X)
				for iParam = 1:optim.NparamOpt
			# Determening if parameters are Log transformed
				if (optim.ParamOpt_LogTransform[iParam]) && !(optim.ParamOpt[iParam]=="Ψm" && optionₘ.σ_2_Ψm == :Constrained)
					Paramₐ = expm1(X[iParam])
				else
					Paramₐ = X[iParam]
				end  # if: optim.ParamOpt_LogTransform

			# Getting the current values of every layer of the hydro parameter of interest
				vectParam = getfield(hydro, Symbol(optim.ParamOpt[iParam]))

			# Updating the value of the parameters for the soil wanting to optimize by keeping the values constant
				vectParam[iZ] = Paramₐ

			# Putting the updated hydro into hydro
				setfield!(hydro, Symbol(optim.ParamOpt[iParam]), vectParam)
		end # for loop

		# ==================== SPECIAL CASE ====================

		# RELATIONSHIP BETWEEN σ AND Ψm
		if (optionₘ.σ_2_Ψm ≠ :No) && ("Ψm" ∈ optim.ParamOpt)
			hydro = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydro, iZ, optionₘ; Pσ=3.0)
		end # optionₘ.σ_2_Ψm ≠ :No

		#  <>=<>=<>=<>=<>=<> Relationship between σ and θr
		if optionₘ.θrOpt==:σ_2_θr && ("θr" ∉ optim.ParamOpt) && ("σ" ∈ optim.ParamOpt)
			hydro.θr[iZ] = hydroRelation.σ_2_θr(hydro, iZ)
		end

		# Converting θsMacMat_ƞ -> θsMacMat
		if  optionₘ.HydroModel == :Kosugi
			hydro.θsMacMat[iZ] = hydro.θsMacMat_ƞ[iZ] * hydro.θs[iZ]
		end

	return hydro
	end  # function: PARAM

end  # module hypixOpt
# ............................................................