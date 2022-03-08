# =============================================================
#		MODULE: residual
# =============================================================
module richard
	import ..timeStep, ..flux, ..ponding, ..residual
	import ..wrc: Ψ_2_θDual, ∂θ∂Ψ, θ_2_ΨDual
	using LinearAlgebra

	export RICHARD_ITERATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD_ITERATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, iCount_ReRun::Int64, discret, Flag_NoConverge::Bool, hydro, iNonConverge::Int64, iT::Int64, IterCount::Int64, NiZ::Int64, paramHypix, Q, Residual, Sorptivity::Float64, Hpond::Vector{Float64}, ΔLnΨmax::Vector{Float64}, ΔPr::Vector{Float64}, ΔSink::Matrix{Float64}, ΔT, Δθₘₐₓ_η::Vector{Float64}, θ::Matrix{Float64}, Ψ::Matrix{Float64}, Ψ_Max::Vector{Float64}, Ψbest::Vector{Float64}, optionHypix)

			# INITIALIZING
			@inbounds @simd for iZ = 1:NiZ
					Ψ[iT,iZ] = Ψ[iT-1,iZ]
				end # for iZ = 1:NiZ
	
			# ITTERATION
			Residual_Max_Best = Inf
			iTer = 0::Int64
			while iTer ≤ paramHypix.N_Iter - 1	
				iTer += 1
				IterCount += 1 # Counting the iterations

				# RESIDUAL MAX BEST: Deriving the Residual max because may be Ψ[iT-1,iZ] is the best solution
				∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, Hpond, θ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, hydro, iT, iTer, NiZ, optionHypix, paramHypix, Q, Residual, Sorptivity, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max)

				# Computing Residual_Max_Best
				if iTer == 1
					Residual_Max_Best = CONVERGENCECRITERIA(discret, iT, NiZ, optionHypix, Residual, ΔT)
				end

				Ψ = SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT, iTer, NiZ, optionHypix, paramHypix, Residual, ΔLnΨmax, Δθₘₐₓ_η, θ, Ψ, Ψ_Max)

				# Averaging the Residuals, depending on method
					Residual_Max = CONVERGENCECRITERIA(discret, iT, NiZ, optionHypix, Residual, ΔT)

				# Determine if iteration made improvement
					if Residual_Max < Residual_Max_Best	
						@inbounds @simd for iZ=1:NiZ
							Ψbest[iZ] = Ψ[iT,iZ]
						end
						Residual_Max_Best = Residual_Max
					end # Residual_Max < Residual_Max_Best 	

				# Did we achieve the goals
				if Residual_Max ≤ paramHypix.WaterBalanceResidual_Max
					break # Move out the loop
				end  # if: Residual
			end # while: iTer ======================

			# Making sure we get the best if convergence fails
			# No convergence
			if iTer == paramHypix.N_Iter
				Flag_NoConverge = true

				# Put the best values
				@inbounds @simd for iZ=1:NiZ
					Ψ[iT,iZ] = Ψbest[iZ]
				end
			else
				Flag_NoConverge = false
			end #  iTer == paramHypix.N_Iter

			for iZ=1:NiZ
				θ[iT,iZ] = Ψ_2_θDual(optionHypix, Ψ[iT,iZ], iZ, hydro)
			end

			# Determine if the simulation is going to rerun with a different time step
			Flag_ReRun, iCount_ReRun, iNonConverge, ΔT = RERUN_HYPIX(discret, Flag_NoConverge, Hpond, hydro, iCount_ReRun, iNonConverge, iT, NiZ, optionHypix, paramHypix, Q, ΔLnΨmax, ΔPr, ΔSink, ΔT, Δθₘₐₓ_η, θ, Ψ)

		return iCount_ReRun, Flag_NoConverge, Flag_ReRun, iNonConverge, iTer,IterCount, Q, Hpond, ΔT, θ, Ψ
		end  # function: RICHARD_SOLVING
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge::Bool, hydro, iT::Int64, iTer::Int64, NiZ::Int64, optionHypix, paramHypix, Q, Residual, Sorptivity, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max)

			if optionHypix.Ponding
				Hpond = ponding.PONDING_SORPTIVITY!(discret, hydro, iT, optionHypix, paramHypix, Q, Sorptivity, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ)
			else
				Hpond[iT] = 0.0::Float64
			end

	#----------------------------------------------------------------
			# ∂R∂Ψ2 = fill(0.0, NiZ)
			# ∂R∂Ψ▽2 = fill(0.0, NiZ)
			# ∂R∂Ψ△2 =  fill(0.0, NiZ)

			for iZ=1:NiZ
				Q, Residual, θ = residual.RESIDUAL(discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, Q, Residual, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ)

				if optionHypix.∂R∂Ψ_NumericalAuto
					∂R∂Ψ[iZ] = residual.∂R∂Ψ_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,NiZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,NiZ)], Ψ_Max)

					∂R∂Ψ▽[iZ]  = residual.∂R∂Ψ▽_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,NiZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,NiZ)], Ψ_Max)
		
					∂R∂Ψ△[iZ]  = residual.∂R∂Ψ△_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,NiZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,NiZ)], Ψ_Max)
				else
					∂R∂Ψ[iZ], ∂R∂Ψ△[iZ], ∂R∂Ψ▽[iZ] = residual.∂RESIDUAL∂Ψ(∂K∂Ψ, discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, ΔT, θ, Ψ)
				end # if optionHypix.∂R∂Ψ_NumericalAuto"
			end #for iZ= 1:NiZ

			# # FOR TESTING...
				# println("One:=================")
				# println("∂R∂Ψ_Deriv=" , ∂R∂Ψ[1:NiZ],"\n") # No good at cell N
				# println("∂R∂Ψ_Num=" , ∂R∂Ψ2[1:NiZ],"\n")

				# println("Two: =================")
				# println("∂R∂Ψ▽_Num=" , ∂R∂Ψ▽[1:NiZ],"\n")
				# println("∂R∂Ψ▽_Der=" , ∂R∂Ψ▽2[1:NiZ],"\n") # No good

				# println("Tree: =================")
				# println("∂R∂Ψ△_Num=" , ∂R∂Ψ△[1:NiZ],"\n") # Good
				# println("∂R∂Ψ△_Der=" , ∂R∂Ψ△2[1:NiZ],"\n")

		return ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, Hpond, θ
		end # function RICHARD
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERGENCECRITERIA
	#     Averaging the Residuals, depending on method
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CONVERGENCECRITERIA(discret, iT::Int64, NiZ::Int64, optionHypix, Residual, ΔT)
			Residual_Norm = 0.0::Float64

			# Does not take into consideration the last cell which has a perfect mass balance
			for iZ = 1:NiZ
				Residual_Norm += (Residual[iZ] / (ΔT[iT] * discret.ΔZ[iZ])) ^ 2.0
			end # for: iZ=NiZ

		return  √(Residual_Norm / Float64(NiZ))			
		end  # function: CONVERGENCECRITERIA
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SOLVING_TRIAGONAL_MATRIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT::Int64, iTer::Int64, NiZ::Int64, optionHypix, paramHypix, Residual, ΔLnΨmax, Δθₘₐₓ_η, θ, Ψ, Ψ_Max)

			Matrix_Trid = Tridiagonal(∂R∂Ψ△[2:NiZ], ∂R∂Ψ[1:NiZ], ∂R∂Ψ▽[1:NiZ-1])

			Residual = reshape(Residual, NiZ, 1) # Transforming from row to column

			NewtonStep = Matrix_Trid \ -Residual
			for iZ=1:NiZ
				# Iteration k-1
					Ψ₀ = Ψ[iT,iZ]
					θ₀ = θ[iT,iZ]
				
				# Updating Ψ
				if !isnan(NewtonStep[iZ])
					# Newton step
					Ψ[iT,iZ] += NewtonStep[iZ]

					# Assuring that the limits of Ψ are physical
						paramHypix.Ψ_MinMin = -10E5
						Ψ[iT,iZ] = min(max(Ψ[iT,iZ], paramHypix.Ψ_MinMin), Ψ_Max[iZ])

					# Correction of θ entering a dry soil 
						if optionHypix.ZhaWetingDrySoil
							Ψ = ZHA_WETING_DRYSOIL(hydro, iT, iZ, optionHypix, θ, θ₀, Ψ, Ψ₀)
						end

					# Smootening the steps
						Δθₘₐₓ_η, Ω = DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT, iZ, optionHypix, paramHypix, ΔLnΨmax, Δθₘₐₓ_η, θ₀, Ψ)
					
					# if optionHypix.DynamicNewtonRaphsonStep
						Ψ[iT,iZ] = Ω * Ψ[iT,iZ] + (1.0 - Ω) * Ψ₀
					# else
					# 	Ψ[iT,iZ] = 0.5 * (Ψ[iT,iZ] +  Ψ₀)

					# end # if optionHypix.DynamicNewtonRaphsonStep

					# if optionHypix.IterReduceOverShoting
					# 	Ψ = Ψ_REDUCE_OVERSHOOTING(iT, iZ, ΔLnΨmax, Ψ, Ψ₀)
					# end
				else
					Ψ[iT,iZ] = Ψ₀ + eps(100.0)
	
				end
			end # for iZ=1:NiZ	

		return Ψ
		end  # function: SOLVING_TRIAGONAL_MATRIX
	# ------------------------------------------------------------------

	
	# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# #		FUNCTION : Ψ_REDUCE_OVERSHOOTING
	# # 		Making sure that the steps of NR are not too big and within the limits of ΔLnΨmax
	# # 		Does not work
	# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# 	function Ψ_REDUCE_OVERSHOOTING(iT, iZ, ΔLnΨmax, Ψ, Ψ₀)
	# 		if Ψ[iT,iZ] ≤ Ψ₀
	# 			Ψ[iT,iZ] = max(Ψ[iT,iZ], expm1(log1p(Ψ₀) - ΔLnΨmax[iZ]))
	# 		else
	# 			Ψ[iT,iZ] = min(Ψ[iT,iZ], expm1(log1p(Ψ₀) + ΔLnΨmax[iZ]))
	# 		end	
	# 	return Ψ
	# 	end  # function: Ψ_REDUCE_OVERSHOOTING
	# #---------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NEWTO_NRAPHSON_STEP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT::Int64, iZ::Int64, optionHypix, paramHypix, ΔLnΨmax, Δθₘₐₓ_η, θ₀, Ψ)

			θ₁ = Ψ_2_θDual(optionHypix, Ψ[iT,iZ], iZ, hydro)

			Δθ = abs(θ₁ - θ₀)

			Δθₘₐₓ = timeStep.ΔθMAX(hydro, iT, iZ, optionHypix, ΔLnΨmax, Ψ) 

			Δθₘₐₓ_η[iZ] = Δθ / Δθₘₐₓ

			Ω =  paramHypix.NewtonStep_Max - (paramHypix.NewtonStep_Max - paramHypix.NewtonStep_Min) * min(Δθₘₐₓ_η[iZ], 1.0) ^ paramHypix.NewtonStep_Power
			
		return Δθₘₐₓ_η, Ω 
		end  # function: NEWTO_NRAPHSON_STEP
	# ------------------------------------------------------------------
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OVERSHOTTING_WET_DRY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	"""Zha, Y., Yang, J., Yin, L., Zhang, Y., Zeng, W., Shi, L., 2017. A modified Picard iteration scheme for overcoming numerical difficulties of simulating infiltration into dry soil. Journal of Hydrology 551, 56–69. https://doi.org/10.1016/j.jhydrol.2017.05.053 """
			function ZHA_WETING_DRYSOIL(hydro, iT, iZ, optionHypix, θ, θ₀, Ψ, Ψ₀)
				# Ψwet = max( 3.5391 * hydro.σ[iZ]^3 - 20.676 * hydro.σ[iZ]^2 + 24.835 * hydro.σ[iZ] + 15.976, 0.0 )

				Ψwet = max(-2.3116 * hydro.σ[iZ] ^ 2.0 - 2.9372 * hydro.σ[iZ] + 27.83, 0.0)

				Ψdry = exp(1.6216 * log(hydro.σ[iZ]) + 8.7268)

				# Determine if there is any oscilation at the wet or dry end of the θ(Ψ) curve
				if (Ψ[iT,iZ] ≤ Ψwet && Ψ₀ ≥ Ψdry)
					θ[iT,iZ] = θ₀ + (Ψ[iT,iZ] - Ψ₀) * ∂θ∂Ψ(optionHypix, Ψ₀, iZ, hydro)

					θ[iT,iZ] = max(min(θ[iT,iZ], hydro.θs[iZ]), hydro.θr[iZ])

					Ψ[iT,iZ] = θ_2_ΨDual(optionHypix, θ[iT,iZ] , iZ, hydro)
				end  # Ψ[iT,iZ] ≤ Ψwet && Ψ₀ ≥ Ψdry
			return Ψ
			end  # function:ZHA_WETING_DRYSOIL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RERUN_HYPIX
	# 		WITH UPDATED Ψ
	#     Rerun if updated ΔT is smaller compared to previously Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RERUN_HYPIX(discret, Flag_NoConverge::Bool, Hpond::Vector{Float64}, hydro, iCount_ReRun::Int64, iNonConverge::Int64, iT::Int64, NiZ::Int64, optionHypix, paramHypix, Q::Matrix{Float64}, ΔLnΨmax::Vector{Float64}, ΔPr::Vector{Float64}, ΔSink::Matrix{Float64}, ΔT::Vector{Float64}, Δθₘₐₓ_η::Vector{Float64}, θ::Matrix{Float64}, Ψ::Matrix{Float64}; Nrerun=2 )

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function RERUN_TIMESTEP⍰(hydro, iT::Int64, NiZ::Int64, optionHypix, ΔLnΨmax::Vector{Float64}, Δθₘₐₓ_η::Vector{Float64}, θ::Matrix{Float64}, Ψ::Matrix{Float64})

					Δθ_Δθₘₐₓ = 0.0::Float64

					for iZ=1:NiZ
						Δθₘₐₓ = timeStep.ΔθMAX(hydro, iT, iZ, optionHypix, ΔLnΨmax, Ψ)
						Δθₘₐₓ_η[iZ] = abs(θ[iT, iZ] -  θ[iT-1, iZ])  / Δθₘₐₓ

						Δθ_Δθₘₐₓ += Δθₘₐₓ_η[iZ]
					end
					Δθ_Δθₘₐₓ = Δθ_Δθₘₐₓ / Float64(NiZ)
			
					if Δθ_Δθₘₐₓ > 1.1 # paramHypix.ΔT_Rerun # paramHypix.ΔT_Rerun
						return true
					else
						return false
					end
				end  # function: RERUN_STEP
			# -------------------------------------------------

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function COMPUTE_ΔT( discret, hydro, iT::Int64, NiZ::Int64, optionHypix, paramHypix, Q::Matrix{Float64}, Hpond::Vector{Float64}, ΔLnΨmax::Vector{Float64}, ΔPr::Vector{Float64}, ΔSink::Matrix{Float64}, ΔT::Vector{Float64}, θ::Matrix{Float64}, Ψ::Matrix{Float64})

					Q[iT,1] = flux.Q!(optionHypix, discret, hydro, 1, iT, NiZ, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT,1], Ψ[iT,1])
					for iZ=1:NiZ
						Q[iT,iZ+1] = flux.Q!(optionHypix, discret, hydro, iZ+1, iT, NiZ, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, min(iZ+1, NiZ)], Ψ[iT,iZ])
					end

					ΔT_New, ~ = timeStep.ADAPTIVE_TIMESTEP(discret, hydro, iT, NiZ, optionHypix, paramHypix, Q, ΔLnΨmax, ΔSink, θ, Ψ)
				return ΔT_New
				end  # function: COMPUTE_ΔT  
			# ------------------------------------------------------------------


			if optionHypix.Flag_ReRun && iCount_ReRun ≤ Nrerun	

				Flag_ReRun = RERUN_TIMESTEP⍰(hydro, iT, NiZ, optionHypix,  ΔLnΨmax, Δθₘₐₓ_η, θ, Ψ)

				ΔTₒ = COMPUTE_ΔT(discret, hydro, iT, NiZ, optionHypix, paramHypix, Q, Hpond, ΔLnΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

				if ΔTₒ * 1.5 < ΔT[iT]
					Flag_ReRun = true
				else
					Flag_ReRun = false
				end

				if Flag_ReRun
					ΔTₒ = COMPUTE_ΔT(discret, hydro, iT, NiZ, optionHypix, paramHypix, Q, Hpond, ΔLnΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

					if ΔTₒ < ΔT[iT]
						ΔT[iT] = ΔTₒ
						iCount_ReRun += 1
						Flag_ReRun = true
					else
						Flag_ReRun = false
					end				
				
				elseif Flag_NoConverge
					Flag_ReRun = true
					iCount_ReRun += 1
					ΔT[iT] =  (ΔT[iT] +  paramHypix.ΔT_Min) * 0.7

				else # <>=<>=<>=<>=<>
					Flag_ReRun = false
					iCount_ReRun = 0
				end
			else
				Flag_ReRun = false
				iCount_ReRun = 0

				if Flag_NoConverge
					iNonConverge += 1
				end
			end  # if: paramHypix.ΔT_Rerun

		return Flag_ReRun, iCount_ReRun, iNonConverge, ΔT
		end  # function: RERUN_HYPIX
	# ------------------------------------------------------------------

		
end # module: richard
#......................................................................