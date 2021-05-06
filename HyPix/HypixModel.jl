# =============================================================
#		module: hypix
# =============================================================
module hypixModel

	import ..checkerror, ..timeStep, ..discretization, ..evapo, ..flux, ..interception, ..interpolate, ..kunsat, ..memory, ..ofHypix, ..option, ..param, ..path, ..pet, ..plot, ..pond, ..priorprocess, ..residual, ..richard, ..rootwateruptake, ..sorptivity, ..tool, ..wrc, ..Δtchange, ..ΨminΨmax

	export HYPIX

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYPIX(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, hydro, Laiᵀ, Laiᵀ_η, N_∑T_Climate, N_iZ::Int64, Q, Residual, veg, Z, ΔEvaporation, ΔHpond, ΔΨmax, ΔPet, ΔPr, ΔSink, ΔT, θ, θ_Ini, Ψ, Ψ_Max, Ψ_Min, Ψbest)

		# Vegetation parameters which vary with time
			for iT = 1:clim.N_Climate
				if option.hypix.LookupTable_Lai
					Laiᵀ[iT]  = (veg.Lai_Max - veg.Lai_Min) * Laiᵀ_η[iT] + veg.Lai_Min
				else
					Laiᵀ[iT] = veg.Lai
				end
				if option.hypix.LookUpTable_CropCoeficient
					CropCoeficientᵀ[iT]  = (veg.CropCoeficient_Max - veg.CropCoeficient_Min) * CropCoeficientᵀ_η[iT]  + veg.CropCoeficient_Min
				else
					CropCoeficientᵀ[iT]  = veg.CropCoeficient
				end
			end # for

		if option.hypix.RainfallInterception
			∑Pet_Climate, ∑Pr_Climate, clim = interception.RAINFALL_INTERCEPTION_START(∑Pet_Climate, ∑Pr_Climate, clim, Laiᵀ, veg)
		end
		
		if option.hypix.RootWaterUptake
			# Last cell of rootzone
			N_iRoot = rootwateruptake.rootDistribution.N_IROOT(N_iZ, veg, Z)

			ΔRootDensity = rootwateruptake.rootDistribution.ROOT_DENSITY(discret, N_iRoot, veg, Z)
		else
			ΔRootDensity = 0.0
			N_iRoot = 1
		end # option.hypix.RootWaterUptake

		# if option.hypix.Evaporation 
		# 	N_iEvapo = evapo.N_IEVAPO(N_iZ, veg, Z) # Depth where evaporation can occure
		# end # option.hypix.Evaporation

		# Minimum or maximum Ψ values this is such that ∂θ∂Ψ ≠ 0 which influences the Newton-Raphson method to be removed
			for iZ=1:N_iZ
				Ψ_Max[iZ], Ψ_Min[iZ] = ΨminΨmax.ΨMINΨMAX(hydro.θs[iZ],  hydro.θsMacMat[iZ],  hydro.σ[iZ],  hydro.σMac[iZ], hydro.Ψm[iZ],  hydro.ΨmMac[iZ])
			end  # for iZ=1:N_iZ

		# IF AdaptiveTimeStep
			if option.hypix.AdaptiveTimeStep == :ΔΨ
				ΔΨmax = timeStep.ΔΨMAX(ΔΨmax, hydro, N_iZ)
			end #  option.hypix.AdaptiveTimeStep == :ΔΨ

		# Water balance of residuals computed from equation
			# WaterBalanceResidual_Max = param.hypix.ΔPr_Error * ∑Pr_Climate[clim.N_Climate] / (Z[N_iZ] * ∑T_Climate[clim.N_Climate])
		
		# FIRST TIME STEP
				# INITIALIZE 

         Flag_NoConverge        = false::Bool
         Flag_ReRun             = false::Bool
         Iter_CountTotal        = 0::Int64
         iNonConverge           = 0::Int64
         iT                     = 1::Int64
         iT_CropCoeficient      = 2::Int64
         iT_Lai                 = 2::Int64
         iT_Pet                 = 2::Int64
         iT_Pr                  = 2::Int64
         ΔEvaporation[1]        = 0.0::Float64
         ΔHpond[1]              = 0.0::Float64
         ΔPet[1]                = 0.0::Float64
         ΔPr[1]                 = 0.0::Float64
         ΔSink[1,1:N_iZ]       .= 0.0::Float64
         ΔT[1]                  = 0.0::Float64
         ∑Pet[1]                = 0.0::Float64
         ∑Pr[1]                 = 0.0::Float64
         ∑T[1]                  = 0.0::Float64
         Count_ReRun            = 0::Int64
			 for iZ = 1:N_iZ
            θ[iT,iZ]   = max( min(hydro.θs[iZ], θ_Ini[iZ]), hydro.θr[iZ]) # Just in case
				Ψ[iT,iZ]   = wrc.θ_2_ΨDual(θ_Ini[iZ], iZ, hydro)
				Ψbest[iZ] = Ψ[iT,iZ] 
				Q[iT,N_iZ] = 0.0
			end

			Q[iT,N_iZ+1] = 0.0

		# =+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+
		while true # this controles the time loop

			# INCREASING OR DECREASING THE TIME STEP
				∑T, FlagContinueLoop, iT, ΔT, Δθ_Max = timeStep.TIMESTEP(∑T, discret, Flag_ReRun, hydro, iT, Float64(N_∑T_Climate), N_iZ, Q, ΔΨmax, ΔSink, ΔT, θ, Ψ)

				if FlagContinueLoop == false
					iT = iT - 1
					break # End of simulation
				end

			# DERIVING FORCING DATA ΔPr & ΔPet:
				∑Pr[iT], ΔPr[iT], iT_Pr = interpolate.∑_2_Δ(∑Pr[iT-1], ∑Pr_Climate, ∑T, ∑T_Climate, iT_Pr, clim.N_Climate, Flag_ReRun, iT)

			# POTENTIAL EVAPOTRANSPIRATION
				if option.hypix.RootWaterUptake || option.hypix.Evaporation
					∑Pet[iT], ΔPet[iT], iT_Pet = interpolate.∑_2_Δ(∑Pet[iT-1], ∑Pet_Climate, ∑T, ∑T_Climate, iT_Pet, clim.N_Climate, Flag_ReRun, iT)
				end # option.hypix.RootWaterUptake || option.hypix.Evaporation

				if option.hypix.Evaporation						
					ΔPet_Evap, ΔPet_Transp = pet.BEER_LAMBERT_LAW(iT, Laiᵀ[iT_Pr-1], ΔPet, veg)
				else
					ΔPet_Transp = ΔPet[iT]
					ΔPet_Evap = 0.0
				end
				
			# ROOT WATER UPTAKE MODEL
				if option.hypix.RootWaterUptake
					ΔSink = rootwateruptake.ROOT_WATER_UPTAKE( CropCoeficientᵀ[iT_Pr-1], iT, N_iRoot, veg, ΔPet_Transp, ΔRootDensity, ΔSink, Ψ)					
				end # option.hypix.RootWaterUptake

			# EVAPORATION FROM THE SURFACE WITH HIGHEST Se
				if option.hypix.Evaporation
					ΔEvaporation = evapo.EVAPORATION!(hydro, iT, ΔEvaporation, ΔPet_Evap, θ)
					
					ΔSink[iT,1] += ΔEvaporation[iT]
				end # option.hypix.Evaporation

			# Checking that not too much water is removed from the layer
				if option.hypix.RootWaterUptake || option.hypix.Evaporation
					for iZ=1:N_iRoot
						ΔSink[iT,iZ] = min(ΔSink[iT,iZ], discret.ΔZ[iZ] * (θ[iT-1,iZ] - hydro.θr[iZ]))
					end
				end # if: option

			# SORPTIVITY TO COMPUTE INFILTRATION RATE				
				Sorptivity = sorptivity.SORPTIVITY(θ[iT-1, 1], 1, hydro; Rtol = 10^-3.0, SorptivityModelScaled = false)

			# SOLVING THE EXPLICIT RICHARDS
				Count_ReRun, Flag_NoConverge, Flag_ReRun, iNonConverge, Iter_CountTotal, Q, ΔHpond, ΔT, θ, Ψ = richard.RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Count_ReRun, discret, Flag_NoConverge, hydro, iNonConverge, iT, Iter_CountTotal, N_iZ, Q, Residual, Sorptivity, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, Δθ_Max, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest)

		end # while loop
		# =+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+	

		N_iT = iT # Maximum of time steps

		return ∑Pet, ∑Pr, ∑T, ∑T_Climate, clim, discret, iNonConverge, Iter_CountTotal, N_iRoot, N_iT, N_iZ, Q, veg, ΔEvaporation, ΔHpond, ΔRootDensity, ΔT, θ, Ψ
	end  # function: HYPIX
	
end  # module hypix
# ............................................................