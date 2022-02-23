# =============================================================
#		module: hypixOpt
# =============================================================
module hypixOpt
	import ..horizonLayer, ..hydroRelation, ..hypixModel, ..ofHypix, ..tool
	using BlackBoxOptim

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIXOPT_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYPIXOPTIMISATION_START(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, hydro, hydro_best, hydroHorizon, hydroHorizon_best, iOpt_Count, Laiᵀ, Laiᵀ_η, Layer, N_∑T_Climate, N_Layer, NiZ, obsTheta, optim, optionHypix, paramHypix, Q, Residual, veg, veg_best, WofBest, Z, ΔEvaporation, Hpond, ΔPet, ΔPr, ΔSink, ΔT, ΔLnΨmax, θ, θini, θSim, Ψ, Ψini, Ψ_Max, Ψ_Min, Ψbest)

		SearchRange = SEARCHRANGE(optim, optionHypix)

		Optimization = BlackBoxOptim.bboptimize(X -> OF_HYPIX(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, hydro, hydroHorizon, Laiᵀ, Laiᵀ_η, Layer, N_∑T_Climate, N_Layer, NiZ, obsTheta, optim, optionHypix, paramHypix, Q, Residual, veg, X, Z, ΔEvaporation, Hpond, ΔPet, ΔPr, ΔSink, ΔT, ΔLnΨmax, θ, θini, θSim, Ψ, Ψini, Ψ_Max, Ψ_Min, Ψbest); SearchRange=SearchRange, NumDimensions=optim.NparamOpt, TraceMode=:silent, MaxFuncEvals=paramHypix.obsTheta.NmaxFuncEvals)

		X = BlackBoxOptim.best_candidate(Optimization)

		hydro, hydroHorizon, veg = PARAM_2_hydro_veg(hydro, hydroHorizon, Layer, N_Layer, NiZ, optim, optionHypix, paramHypix, X, veg)

		WofBest[iOpt_Count] = BlackBoxOptim.best_fitness(Optimization)

		println("\n			~ WOFbest = ", WofBest[iOpt_Count] , "\n")

		if iOpt_Count == 1
         hydro_best        = deepcopy(hydro)
         hydroHorizon_best = deepcopy(hydroHorizon)
			veg_best          = deepcopy(veg)

		elseif iOpt_Count ≥ 2 && WofBest[iOpt_Count] <  WofBest[iOpt_Count-1]
			println("\n		== ~ IMPROVING ~ == \n")			
				hydro_best        = deepcopy(hydro) 
            hydroHorizon_best = deepcopy(hydroHorizon)
				veg_best          = deepcopy(veg)
									
		# Sorry Not improving so we revert
		elseif iOpt_Count ≥ 2 && WofBest[iOpt_Count] ≥ WofBest[iOpt_Count-1]
			println("\n		== ~ NO IMPROVEMENTS ~ == \n")
            hydro        = deepcopy(hydro_best)
            hydroHorizon = deepcopy(hydroHorizon_best)
				veg          = deepcopy(veg_best)

				WofBest[iOpt_Count] = WofBest[iOpt_Count-1]			
		end # if iOpt_Count

	return hydro, hydro_best, hydroHorizon, hydroHorizon_best, veg, veg_best, WofBest
	end  # function: HYPIXOPT_START



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_HYPIX(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, hydro, hydroHorizon, Laiᵀ, Laiᵀ_η, Layer, N_∑T_Climate, N_Layer, NiZ, obsTheta, optim, optionHypix, paramHypix, Q, Residual, veg, X, Z, ΔEvaporation, Hpond, ΔPet, ΔPr, ΔSink, ΔT, ΔLnΨmax, θ, θini, θSim, Ψ, Ψini, Ψ_Max, Ψ_Min, Ψbest)

			# New optimized paramHypix which are put into the matching veg or hydro parameters
				hydro, hydroHorizon, veg = PARAM_2_hydro_veg(hydro, hydroHorizon, Layer, N_Layer, NiZ, optim, optionHypix, paramHypix, X, veg)
		
			# Running Hypix model	
				∑Pet, ∑Pr, ∑T, ∑T_Climate, clim, discret, iNonConverge, IterCount, N_iRoot, Nit, NiZ, Q, veg, ΔEvaporation, Hpond, ΔRootDensity, ΔT, θ, Ψ  = hypixModel.HYPIX(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, hydro, Laiᵀ, Laiᵀ_η, N_∑T_Climate, NiZ, optionHypix, paramHypix, Q, Residual, veg, Z, ΔEvaporation, Hpond, ΔPet, ΔPr, ΔSink, ΔT, ΔLnΨmax, θ, θini, Ψ, Ψini, Ψ_Max, Ψ_Min, Ψbest)

			# Weighted Objective Function
				Wof = ofHypix.WOF_θ(∑T[1:Nit], Nit, NiZ, obsTheta, paramHypix, Hpond[1:Nit], θ[1:Nit,1:NiZ], θSim)
				
		println("\n			~ Wof = $Wof ~")
		return Wof
		end  # function: OF_HYPIX


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SEARCHRANGE
	#		Required by BlackBoxOptim
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SEARCHRANGE(optim, optionHypix)
			ParamOpt_Min₂ = copy(optim.ParamOpt_Min)
			ParamOpt_Max₂ = copy(optim.ParamOpt_Max)

			# Making sure that for constrained optimisation Ψm is between 0 & 1
			if (optionHypix.σ_2_Ψm⍰=="Constrained") && ("Ψm" ∈ optim.ParamOpt)
				iψm = findfirst(isequal("Ψm"), optim.ParamOpt)[1]

				ParamOpt_Min₂[iψm] = 0.0
				ParamOpt_Max₂[iψm] = 1.0
			end # optionHypix.σ_2_Ψm⍰==Constrained

      # "θs_Opt⍰"                 = "No" #  <θs_Opt> θs is derived by multiplying a parameter to Max(θobs) for all profiles; <No>
			if  ("θs" ∈ optim.ParamOpt) && (optionHypix.θs_Opt⍰ ≠ "No")
				iθs = findfirst(isequal("θs"), optim.ParamOpt)[1]

				ParamOpt_Min₂[iθs] = 0.0
				ParamOpt_Max₂[iθs] = 1.0
			end # "θs" ∈ optim.ParamOpt

			SearchRange = (collect(zip(Float64.(ParamOpt_Min₂), Float64.(ParamOpt_Max₂))))
			println("	~ SearchRange = ", SearchRange)
		return SearchRange
		end  # function: SEARCHRANGE



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PARAM_2_hydro_veg
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PARAM_2_hydro_veg(hydro, hydroHorizon, Layer, N_Layer, NiZ, optim, optionHypix, paramHypix, X, veg)

			println("\n		==== OPT PARAM ====")
			for iParam = 1:optim.NparamOpt

				# Determening if parameters are Log transformed
				if optim.ParamOpt_LogTransform[iParam]
					Paramₐ = expm1(X[iParam])
				else
					Paramₐ = X[iParam]
				end  # if: optim.ParamOpt_LogTransform


				if optim.ParamOpt_Type[iParam] == "hydro"
					# Getting the current values of every layer of the hydro parameter of interest
						vectParam = getfield(hydroHorizon, Symbol(optim.ParamOpt[iParam]))
					
					# Horizons wanting to optimize of the selected hydraulic parameter
						iHorizon_Start = optim.ParamOpt_HorizonEq[iParam][1]
						iHorizon_End   = optim.ParamOpt_HorizonEq[iParam][2]
					
					# Updating the value of the parameters for the horizons wanting to optimize by keeping the values constant
						for iZ = iHorizon_Start:iHorizon_End
							vectParam[iZ] = Paramₐ
						end  # for iZ

					# Putting the updated hydro paramHypix into hydrohydroHorizon
						setfield!(hydroHorizon, Symbol(optim.ParamOpt[iParam]), vectParam)

				elseif optim.ParamOpt_Type[iParam] == "veg"
					# Putting veg paramHypix in hydroHorizon
					setfield!(veg, Symbol(optim.ParamOpt[iParam]), Paramₐ)

				end # if optim.ParamOpt_Type[iParam]

			end # for loop

			# ==================== SPECIAL CASE ====================

			# RELATIONSHIP BETWEEN σ AND Ψm
			if (optionHypix.σ_2_Ψm⍰ ≠ "No") && ("Ψm" ∈ optim.ParamOpt)
		
				# <>=<>=<>=<>=<>=<> Horizons wanting to optimize the selected hydraulic parameter
					iParam = findfirst(isequal("σ"), optim.ParamOpt)[1]

					iHorizon_Start = optim.ParamOpt_HorizonEq[iParam][1]
					iHorizon_End   = optim.ParamOpt_HorizonEq[iParam][2]

					hydroHorizon = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydroHorizon, iHorizon_Start, optionHypix, paramHypix)

				# Updating the horizons which are optimised simultaneously
					for iZ = iHorizon_Start:iHorizon_End
						hydroHorizon.Ψm[iZ] = hydroHorizon.Ψm[iHorizon_Start]
					end  # for iZ
			end # optionHypix.σ_2_Ψm⍰ ≠ No

			#  <>=<>=<>=<>=<>=<> Relationship between σ and θr
				if optionHypix.σ_2_θr && ("θr" ∉ optim.ParamOpt) && ("σ" ∈ optim.ParamOpt)
					iParam = findfirst(isequal("σ"), optim.ParamOpt)[1]

					iHorizon_Start = optim.ParamOpt_HorizonEq[iParam][1]
					iHorizon_End   = optim.ParamOpt_HorizonEq[iParam][2]
	
				# Updating the horizons which are optimised simultaneously
					for iZ = iHorizon_Start:iHorizon_End
						hydroHorizon.θr[iZ] = hydroRelation.σ_2_θr(hydroHorizon, iZ)
					end # iZ	
				end

			#  <>=<>=<>=<>=<>=<> Assuring the limits of 
				if  ("θs" ∈ optim.ParamOpt) && (optionHypix.θs_Opt⍰ == "θs_Opt")
					for iZ = iHorizon_Start:iHorizon_End
						hydroHorizon.θs[iZ] = tool.norm.∇NORM_2_PARAMETER(hydroHorizon.θs[iZ], hydroHorizon.θs_Min[iZ], hydroHorizon.θs_Max[iZ])
					end # iZ
				end # if  ("θs" ∈ optim.ParamOpt) && (optionHypix.θs_Opt⍰ == :θs_Opt)


			#  <>=<>=<>=<>=<>=<> Assuring the limits of θs are physical
				if  ("θs" ∈ optim.ParamOpt)
					for iZ = iHorizon_Start:iHorizon_End
						hydroHorizon.θs[iZ] = min( max(hydroHorizon.θs[iZ], hydroHorizon.θs_Min[iZ]), hydroHorizon.θs_Max[iZ])
					end # iZ
				end # if  ("θs" ∈ optim.ParamOpt) 
		

			# Converting θsMacMat_ƞ -> θsMacMat
				for iZ=1:N_Layer
					hydroHorizon.θsMacMat[iZ] = hydroHorizon.θsMacMat_ƞ[iZ] * hydroHorizon.θs[iZ]
				end 

				if "hydro" ∈ optim.ParamOpt_Type
					println("\n			~ θs = ",  hydroHorizon.θs[iHorizon_Start])
					println("			~ θsMacMat = ",  hydroHorizon.θsMacMat[iHorizon_Start])
					println("			~ θr = ",  hydroHorizon.θr[iHorizon_Start])
					println("			~ Ψm = ",  hydroHorizon.Ψm[iHorizon_Start])
					println("			~ σ = ",  hydroHorizon.σ[iHorizon_Start])
					println("			~ Ks = ",  hydroHorizon.Ks[iHorizon_Start])
				end

			# Transforming horizonLayer -> hydro
				hydro = horizonLayer.HYDROHORIZON_2_HYDRO(hydroHorizon, Layer, NiZ, optionHypix)

		return hydro, hydroHorizon, veg
		end  # function: PARAM_2_hydro_veg
		
	end  # module hypixOpt
# ............................................................rm 