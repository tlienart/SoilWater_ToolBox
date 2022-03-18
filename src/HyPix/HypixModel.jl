# =============================================================
#		module: hypix
# =============================================================
module hypixModel

	import ..evaporation, ..interception, ..interpolate, ..pet, ..richard, ..rootWaterUptake, ..sorptivity, ..timeStep, ..ΨminΨmax
	import ..wrc: θ_2_ΨDual, Ψ_2_θDual

	export HYPIX_MODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIX_MODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYPIX_MODEL(∂K∂Ψ::Vector{Float64}, ∂R∂Ψ::Vector{Float64}, ∂R∂Ψ△::Vector{Float64}, ∂R∂Ψ▽::Vector{Float64}, ∑Pet::Vector{Float64}, ∑Pet_Climate::Vector{Float64}, ∑Pr::Vector{Float64}, ∑Pr_Climate::Vector{Float64}, ∑T::Vector{Float64}, ∑T_Climate::Vector{Float64}, clim, CropCoeficientᵀ::Vector{Float64}, CropCoeficientᵀ_η::Vector{Float64}, discret, Flag_θΨini::Symbol, hydro, Laiᵀ::Vector{Float64}, Laiᵀ_η::Vector{Float64}, N_∑T_Climate::Int64, NiZ::Int64, optionHypix, paramHypix, Q::Matrix{Float64}, Residual::Vector{Float64}, veg, Z::Vector{Float64}, ΔEvaporation::Vector{Float64}, Hpond::Vector{Float64}, ΔPet::Vector{Float64}, ΔPr::Vector{Float64}, ΔRunoff::Vector{Float64}, ΔSink::Matrix{Float64}, ΔT::Vector{Float64}, ΔLnΨmax::Vector{Float64}, θ::Matrix{Float64}, θini_or_Ψini::Vector{Float64}, Ψ::Matrix{Float64}, Ψ_Max::Vector{Float64}, Ψ_Min::Vector{Float64}, Ψbest::Vector{Float64})

		# VEGETATION PARAMETERS WHICH VARY WITH TIME
			for iT = 1:clim.N_Climate
				if optionHypix.LookupTable_Lai
					Laiᵀ[iT]  = (veg.Lai_Max - veg.Lai_Min) * Laiᵀ_η[iT] + veg.Lai_Min
				else
					Laiᵀ[iT] = veg.Lai
				end
				if optionHypix.LookUpTable_CropCoeficient
					CropCoeficientᵀ[iT]  = (veg.CropCoeficient_Max - veg.CropCoeficient_Min) * CropCoeficientᵀ_η[iT]  + veg.CropCoeficient_Min
				else
					CropCoeficientᵀ[iT]  = veg.CropCoeficient
				end
			end # for

		# RAINFALL INTERCEPTION
		if optionHypix.RainfallInterception
			∑Pet_Climate, ∑Pr_Climate, clim = interception.RAINFALL_INTERCEPTION_START(∑Pet_Climate, ∑Pr_Climate, clim, Laiᵀ, optionHypix, veg)
		end
		
		# ROOTS
		if optionHypix.RootWaterUptake
			N_iRoot = rootWaterUptake.rootDistribution.N_IROOT(NiZ, veg, Z) # Last cell of rootzone

			ΔRootDensity = rootWaterUptake.rootDistribution.ROOT_DENSITY(discret, N_iRoot, veg, Z)
		else
			ΔRootDensity = 0.0::Float64
			N_iRoot = 1::Int64
		end # optionHypix.RootWaterUptake

		# if optionHypix.Evaporation 
		# 	N_iEvapo = evaporation.N_IEVAPO(NiZ, veg, Z) # Smap_Depth where evaporation can occure
		# end # optionHypix.Evaporation

		# MINIMUM OR MAXIMUM Ψ VALUES THIS IS SUCH THAT ∂Θ∂Ψ ≠ 0 WHICH INFLUENCES THE NEWTON-RAPHSON METHOD TO BE REMOVED
			for iZ=1:NiZ
            Ψ_Max[iZ],~ = ΨminΨmax.ΨMINΨMAX(hydro.θs[iZ], hydro.θsMacMat[iZ],  hydro.σ[iZ],  hydro.σMac[iZ], hydro.Ψm[iZ], hydro.ΨmMac[iZ])
			end  # for iZ=1:NiZ

		# ADAPTIVETIMESTEP
			ΔLnΨmax = timeStep.ΔΨMAX!(hydro, NiZ, optionHypix, paramHypix, ΔLnΨmax)

		# FIRST TIME STEP
         Flag_NoConverge        = false::Bool
         Flag_ReRun             = false::Bool
         IterCount              = 0::Int64
			iTer = 10
         iNonConverge           = 0::Int64
         iT                     = 1::Int64
         iT_Pet                 = 2::Int64
         iT_Pr                  = 2::Int64
         ΔEvaporation[1]        = 0.0::Float64
         Hpond[1]              = 0.0::Float64
         ΔPet[1]                = 0.0::Float64
         ΔPr[1]                 = 0.0::Float64
         ΔSink[1,1:NiZ]        .= 0.0::Float64
         ΔT[1]                  = 0.0::Float64
         ∑Pet[1]                = 0.0::Float64
         ∑Pr[1]                 = 0.0::Float64
         ∑T[1]                  = 0.0::Float64
         iCount_ReRun            = 1::Int64
			
		# Boundary conditions
			if Flag_θΨini == :θini
				for iZ = 1:NiZ
               θ[1,iZ] = max( min(hydro.θs[iZ], θini_or_Ψini[iZ]), hydro.θr[iZ] ) # Just in case
               Ψ[1,iZ] = θ_2_ΨDual(optionHypix, θini_or_Ψini[iZ], iZ, hydro)
				end

			elseif Flag_θΨini == :Ψini
				for iZ = 1:NiZ
               Ψ[1,iZ] = θini_or_Ψini[iZ]
               θ[1,iZ] = Ψ_2_θDual(optionHypix, θini_or_Ψini[iZ], iZ, hydro)
				end
			end

			if optionHypix.TopBoundary⍰ == "Ψ"
				Ψ[1,1] = paramHypix.Ψ_Top
				θ[1,1]  = Ψ_2_θDual(optionHypix, paramHypix.Ψ_Top, 1, hydro)
			end

			if optionHypix.BottomBoundary⍰ == "Ψ"
				Ψ[1,NiZ] = paramHypix.Ψ_Botom
				θ[1,NiZ]  = Ψ_2_θDual(optionHypix, paramHypix.Ψ_Botom, NiZ, hydro)	
			end

			for iZ = 1:NiZ
            Ψbest[iZ] = Ψ[1,iZ]
            Q[1,NiZ]  = 0.0::Float64
			end
			Q[1,NiZ+1] = 0.0::Float64

		# =+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+
		while true # this controles the time loop

			# INCREASING OR DECREASING THE TIME STEP
				∑T, FlagContinueLoop, iT, ΔT, Δθ_Max = timeStep.TIMESTEP(∑T, discret, Flag_ReRun, hydro, iT, iTer, Float64(N_∑T_Climate), NiZ, optionHypix, paramHypix, Q, ΔLnΨmax, ΔSink, ΔT, θ, Ψ)

				if FlagContinueLoop == false
					iT = iT - 1
					break # End of simulation
				end

			# DERIVING FORCING DATA ΔPr & ΔPet:
				∑Pr[iT], ΔPr[iT], iT_Pr = interpolate.∑_2_Δ(∑Pr[iT-1], ∑Pr_Climate, ∑T, ∑T_Climate, iT_Pr, clim.N_Climate, Flag_ReRun, iT)

			# POTENTIAL EVAPOTRANSPIRATION
				if optionHypix.RootWaterUptake || optionHypix.Evaporation
					∑Pet[iT], ΔPet[iT], iT_Pet = interpolate.∑_2_Δ(∑Pet[iT-1], ∑Pet_Climate, ∑T, ∑T_Climate, iT_Pet, clim.N_Climate, Flag_ReRun, iT)
				end # optionHypix.RootWaterUptake || optionHypix.Evaporation

				if optionHypix.Evaporation						
					ΔPet_Evap, ΔPet_Transp = pet.BEER_LAMBERT_LAW(iT, Laiᵀ[iT_Pr-1], ΔPet, veg)
				else
					ΔPet_Transp = ΔPet[iT]
					ΔPet_Evap = 0.0::Float64
				end
				
			# ROOT WATER UPTAKE MODEL
				if optionHypix.RootWaterUptake
					ΔSink = rootWaterUptake.ROOT_WATER_UPTAKE(CropCoeficientᵀ[iT_Pr-1], iT, N_iRoot, optionHypix, veg, ΔPet_Transp, ΔRootDensity, ΔSink, Ψ)					
				end # optionHypix.RootWaterUptake

			# EVAPORATION FROM THE SURFACE WITH HIGHEST Se
				if optionHypix.Evaporation
					ΔEvaporation = evaporation.EVAPORATION!(hydro, iT, ΔEvaporation, ΔPet_Evap, θ)
					
					ΔSink[iT,1] += ΔEvaporation[iT]
				end # optionHypix.Evaporation

			# Checking that not too much water is removed from the layer
				if optionHypix.RootWaterUptake || optionHypix.Evaporation
					for iZ=1:N_iRoot
						ΔSink[iT,iZ] = min(ΔSink[iT,iZ], discret.ΔZ[iZ] * (θ[iT-1,iZ] - hydro.θr[iZ]))
					end
				end # if: optionHypix

			# SORPTIVITY TO COMPUTE INFILTRATION RATE
				Sorptivity = sorptivity.SORPTIVITY(θ[iT-1, 1], 1, hydro, optionHypix, optionHypix; Rtol = 10^-3.0, SorptivityModelScaled=false)
		
			# SOLVING THE EXPLICIT RICHARDS
			Flag_NoConverge, Flag_ReRun, Hpond, iCount_ReRun, iNonConverge, iTer, IterCount, Q, ΔRunoff, ΔT, θ, Ψ = richard.RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, iCount_ReRun, discret, Flag_NoConverge, hydro, iNonConverge, iT, IterCount, NiZ, paramHypix, Q, Residual, Sorptivity, Hpond, ΔLnΨmax, ΔPr, ΔRunoff, ΔSink, ΔT, θ, Ψ, Ψ_Min, Ψ_Max, Ψbest, optionHypix)
				
			# SPECIAL BOUNDARY CONDITIONS
				if optionHypix.TopBoundary⍰ == "Ψ"
					ΔPr[iT] = ΔT[iT] * Q[iT, 1]
					∑Pr[iT] = ∑Pr[iT-1] + ΔPr[iT]
				end
		end # while loop
		# =+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+	

		Nit = iT # Maximum time steps

	return ∑Pet, ∑Pr, ∑T, ∑T_Climate, clim, discret, Hpond, iNonConverge, IterCount, N_iRoot, Nit, NiZ, Q, veg, ΔEvaporation, ΔRootDensity, ΔRunoff, ΔT, θ, Ψ
	end  # function: HYPIX_MODEL
	
end  # module hypix
# ............................................................