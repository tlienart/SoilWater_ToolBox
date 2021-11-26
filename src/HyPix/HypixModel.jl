# =============================================================
#		module: hypix
# =============================================================
module hypixModel

	import ..evaporation, ..interception, ..interpolate, ..pet, ..richard, ..rootWaterUptake, ..sorptivity, ..timeStep, ..ΨminΨmax
	import ..wrc: θ_2_ΨDual, Ψ_2_θDual

	export HYPIX

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYPIX(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, hydro, Laiᵀ, Laiᵀ_η, N_∑T_Climate, NiZ::Int64, option, param, Q, Residual, veg, Z, ΔEvaporation, ΔHpond, ΔPet, ΔPr, ΔSink, ΔT, ΔLnΨmax, θ, θini, Ψ, Ψini, Ψ_Max, Ψ_Min, Ψbest)

		# VEGETATION PARAMETERS WHICH VARY WITH TIME
			for iT = 1:clim.N_Climate
				if option.hyPix.LookupTable_Lai
					Laiᵀ[iT]  = (veg.Lai_Max - veg.Lai_Min) * Laiᵀ_η[iT] + veg.Lai_Min
				else
					Laiᵀ[iT] = veg.Lai
				end
				if option.hyPix.LookUpTable_CropCoeficient
					CropCoeficientᵀ[iT]  = (veg.CropCoeficient_Max - veg.CropCoeficient_Min) * CropCoeficientᵀ_η[iT]  + veg.CropCoeficient_Min
				else
					CropCoeficientᵀ[iT]  = veg.CropCoeficient
				end
			end # for

		# RAINFALL INTERCEPTION
		if option.hyPix.RainfallInterception
			∑Pet_Climate, ∑Pr_Climate, clim = interception.RAINFALL_INTERCEPTION_START(∑Pet_Climate, ∑Pr_Climate, clim, Laiᵀ, option, veg)
		end
		
		# ROOTS
		if option.hyPix.RootWaterUptake
			N_iRoot = rootWaterUptake.rootDistribution.N_IROOT(NiZ, veg, Z)# Last cell of rootzone

			ΔRootDensity = rootWaterUptake.rootDistribution.ROOT_DENSITY(discret, N_iRoot, veg, Z)
		else
			ΔRootDensity = 0.0
			N_iRoot = 1
		end # option.hyPix.RootWaterUptake

		# if option.hyPix.Evaporation 
		# 	N_iEvapo = evaporation.N_IEVAPO(NiZ, veg, Z) # Smap_Depth where evaporation can occure
		# end # option.hyPix.Evaporation

		# MINIMUM OR MAXIMUM Ψ VALUES THIS IS SUCH THAT ∂Θ∂Ψ ≠ 0 WHICH INFLUENCES THE NEWTON-RAPHSON METHOD TO BE REMOVED
			for iZ=1:NiZ
				Ψ_Max[iZ],~ = ΨminΨmax.ΨMINΨMAX(hydro.θs[iZ],  hydro.θsMacMat[iZ],  hydro.σ[iZ],  hydro.σMac[iZ], hydro.Ψm[iZ], hydro.ΨmMac[iZ])
				Ψ_Min[iZ] = param.hyPix.Ψ_MinMin
			end  # for iZ=1:NiZ

		# ADAPTIVETIMESTEP
			ΔLnΨmax = timeStep.ΔΨMAX!(hydro, NiZ, option, param, ΔLnΨmax)

		# FIRST TIME STEP
         Flag_NoConverge        = false::Bool
         Flag_ReRun             = false::Bool
         IterCount              = 0::Int64
         iNonConverge           = 0::Int64
         iT                     = 1::Int64
         iT_Pet                 = 2::Int64
         iT_Pr                  = 2::Int64
         ΔEvaporation[1]        = 0.0::Float64
         ΔHpond[1]              = 0.0::Float64
         ΔPet[1]                = 0.0::Float64
         ΔPr[1]                 = 0.0::Float64
         ΔSink[1,1:NiZ]        .= 0.0::Float64
         ΔT[1]                  = 0.0::Float64
         ∑Pet[1]                = 0.0::Float64
         ∑Pr[1]                 = 0.0::Float64
         ∑T[1]                  = 0.0::Float64
         Count_ReRun            = 0::Int64
			
		# Boundary conditions
			if Flag_θΨini == :θini
				Ψini =  fill(0.0::Float64, NiZ)
				
				for iZ = 1:NiZ
					θ[1,iZ]   = max( min(hydro.θs[iZ], θini[iZ]), hydro.θr[iZ] ) # Just in case
					Ψ[1,iZ]   = θ_2_ΨDual(option.hydro, θini[iZ], iZ, hydro)

					if iZ == 1 && option.hyPix.TopBoundary⍰ == "Ψ"
						Ψ[1,1] = param.hyPix.Ψ_Top
						θ[1,1]  = Ψ_2_θDual(option.hydro, Ψ[iT,1], iZ, hydro)
					end
	
					if iZ == NiZ && option.hyPix.BottomBoundary⍰ == "Ψ"
						Ψ[1,NiZ] = param.hyPix.Ψ_Botom
						θ[1,NiZ]  = Ψ_2_θDual(option.hydro, Ψ[1,NiZ], NiZ, hydro)
					end

					Ψini[iZ] = Ψ[1,iZ]

					Ψbest[iZ]  = Ψini[iZ]
					Q[1,NiZ] = 0.0
				end

			elseif Flag_θΨini == :Ψini
				θini = fill(0.0::Float64, NiZ)

				if option.hyPix.TopBoundary⍰ == "Ψ"
					Ψini[1] = param.hyPix.Ψ_Top
				end

				if option.hyPix.BottomBoundary⍰ == "Ψ"
					Ψini[NiZ] = param.hyPix.Ψ_Botom
				end

				for iZ = 1:NiZ
					Ψ[1,iZ] = Ψini[iZ]
					θ[1,iZ]  = Ψ_2_θDual(option.hydro, Ψini[iZ], iZ, hydro)
					θini[iZ] = θ[1,iZ]

					Ψbest[iZ]  = Ψ[1,iZ]
					Q[1,NiZ] = 0.0
				end
			end

			Q[1,NiZ+1] = 0.0

		# =+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+
		while true # this controles the time loop

			# INCREASING OR DECREASING THE TIME STEP
				∑T, FlagContinueLoop, iT, ΔT, Δθ_Max = timeStep.TIMESTEP(∑T, discret, Flag_ReRun, hydro, iT, Float64(N_∑T_Climate), NiZ, option, param, Q, ΔLnΨmax, ΔSink, ΔT, θ, Ψ)

				if FlagContinueLoop == false
					iT = iT - 1
					break # End of simulation
				end

			# DERIVING FORCING DATA ΔPr & ΔPet:
				∑Pr[iT], ΔPr[iT], iT_Pr = interpolate.∑_2_Δ(∑Pr[iT-1], ∑Pr_Climate, ∑T, ∑T_Climate, iT_Pr, clim.N_Climate, Flag_ReRun, iT)

			# POTENTIAL EVAPOTRANSPIRATION
				if option.hyPix.RootWaterUptake || option.hyPix.Evaporation
					∑Pet[iT], ΔPet[iT], iT_Pet = interpolate.∑_2_Δ(∑Pet[iT-1], ∑Pet_Climate, ∑T, ∑T_Climate, iT_Pet, clim.N_Climate, Flag_ReRun, iT)
				end # option.hyPix.RootWaterUptake || option.hyPix.Evaporation

				if option.hyPix.Evaporation						
					ΔPet_Evap, ΔPet_Transp = pet.BEER_LAMBERT_LAW(iT, Laiᵀ[iT_Pr-1], ΔPet, veg)
				else
					ΔPet_Transp = ΔPet[iT]
					ΔPet_Evap = 0.0
				end
				
			# ROOT WATER UPTAKE MODEL
				if option.hyPix.RootWaterUptake
					ΔSink = rootWaterUptake.ROOT_WATER_UPTAKE(CropCoeficientᵀ[iT_Pr-1], iT, N_iRoot, option, veg, ΔPet_Transp, ΔRootDensity, ΔSink, Ψ)					
				end # option.hyPix.RootWaterUptake

			# EVAPORATION FROM THE SURFACE WITH HIGHEST Se
				if option.hyPix.Evaporation
					ΔEvaporation = evaporation.EVAPORATION!(hydro, iT, ΔEvaporation, ΔPet_Evap, θ)
					
					ΔSink[iT,1] += ΔEvaporation[iT]
				end # option.hyPix.Evaporation

			# Checking that not too much water is removed from the layer
				if option.hyPix.RootWaterUptake || option.hyPix.Evaporation
					for iZ=1:N_iRoot
						ΔSink[iT,iZ] = min(ΔSink[iT,iZ], discret.ΔZ[iZ] * (θ[iT-1,iZ] - hydro.θr[iZ]))
					end
				end # if: option

			# SORPTIVITY TO COMPUTE INFILTRATION RATE				
				Sorptivity = sorptivity.SORPTIVITY(θ[iT-1, 1], 1, hydro, option, option.hydro; Rtol = 10^-3.0, SorptivityModelScaled=false)
		
			# SOLVING THE EXPLICIT RICHARDS
				Count_ReRun, Flag_NoConverge, Flag_ReRun, iNonConverge, IterCount, Q, ΔHpond, ΔT, θ, Ψ = richard.RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Count_ReRun, discret, Flag_NoConverge, hydro, iNonConverge, iT, IterCount, NiZ, param, Q, Residual, Sorptivity, ΔHpond, ΔLnΨmax, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψbest, option)
				
			# SPECIAL BOUNDARY CONDITIONS
				if option.hyPix.TopBoundary⍰ == "Ψ"
					ΔPr[iT] = ΔT[iT] * Q[iT, 1]
					∑Pr[iT] = ∑Pr[iT-1] + ΔPr[iT]
				end
		end # while loop
		# =+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+	

		Nit = iT # Maximum time steps

	return ∑Pet, ∑Pr, ∑T, ∑T_Climate, clim, discret, iNonConverge, IterCount, N_iRoot, Nit, NiZ, Q, veg, ΔEvaporation, ΔHpond, ΔRootDensity, ΔT, θ, Ψ
	end  # function: HYPIX
	
end  # module hypix
# ............................................................