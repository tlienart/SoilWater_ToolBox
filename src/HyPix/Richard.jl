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
		function RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Count_ReRun::Int64, discret, Flag_NoConverge::Bool, hydro, iNonConverge::Int64, iT::Int64, IterCount::Int64, NiZ::Int64, param, Q, Residual, Sorptivity::Float64, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest, option)

			# INITIALIZING
			@inbounds @simd for iZ = 1:NiZ
					Ψ[iT,iZ] = Ψbest[iZ]
				end # for iZ = 1:NiZ

				# Residual_Max_Best it should be improved compared to Ψ[iT,iZ] = Ψ[iT-1,iZ]
				~, ~, ~, ~, Residual, ~, ~ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, hydro, iT, NiZ, param, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, option)

				# Averaging the Residuals, depending on method
					Residual_Max_Best = RESIDUAL_MAX(discret, iT, NiZ, option, Residual, ΔT)
	
			# ITTERATION
			iTer = 0::Int64
			while iTer ≤ param.hyPix.N_Iter	
				iTer += 1
				IterCount += 1 # Counting the iterations
		
				∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, ΔHpond, θ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, hydro, iT, NiZ, param, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, option)

				Ψ = SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT, iTer, NiZ, option, param, Residual, ΔΨmax, θ, Ψ)

				# Averaging the Residuals, depending on method
					Residual_Max = RESIDUAL_MAX(discret, iT, NiZ, option, Residual, ΔT)

				# Determine if iteration made improvement
					if Residual_Max < Residual_Max_Best	
						for iZ=1:NiZ
							Ψbest[iZ] = Ψ[iT,iZ]
						end
						Residual_Max_Best = Residual_Max
					end # Residual_Max < Residual_Max_Best 	

				# Did we achieve the goals
				if Residual_Max ≤ param.hyPix.WaterBalanceResidual_Max
					break # Move out the loop
				end  # if: Residual
			end # while: iTer ======================

			# Making sure we get the best if convergence fails
			if iTer ≥ param.hyPix.N_Iter
				Flag_NoConverge = true

				iNonConverge += 1

				# Put the best values
				for iZ=1:NiZ
					Ψ[iT,iZ] = Ψbest[iZ]
				end
			else
				Flag_NoConverge = false

			end #  iTer == param.hyPix.N_Iter

			# Determine if the simulation is going to rerun with a different time step
				Count_ReRun, Flag_ReRun, ΔT = RERUN_HYPIX(Count_ReRun, discret, Flag_NoConverge, hydro, iT, NiZ, option, param, Q, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

		return Count_ReRun, Flag_NoConverge, Flag_ReRun, iNonConverge, IterCount, Q, ΔHpond, ΔT, θ, Ψ
		end  # function: RICHARD_SOLVING
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge::Bool, hydro, iT::Int64, NiZ::Int64, param, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, option)

			ΔHpond = ponding.PONDING_SORPTIVITY!(discret, hydro, iT, option, param, Q, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ)

	#----------------------------------------------------------------
			# ∂R∂Ψ2 = fill(0.0, NiZ)
			# ∂R∂Ψ▽2 = fill(0.0, NiZ)
			# ∂R∂Ψ△2 =  fill(0.0, NiZ)

			for iZ=1:NiZ		
				Q, Residual, θ = residual.RESIDUAL(option, discret, hydro, iT, iZ, NiZ, param, Q, Residual, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ)

				if option.hyPix.∂R∂Ψ_Numerical
					∂R∂Ψ[iZ] = residual.∂R∂Ψ_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, option, param, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,NiZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,NiZ)], Ψ_Max)

					∂R∂Ψ▽[iZ]  = residual.∂R∂Ψ▽_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, option, param, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,NiZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,NiZ)], Ψ_Max)
		
					∂R∂Ψ△[iZ]  = residual.∂R∂Ψ△_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, option, param, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,NiZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,NiZ)], Ψ_Max)
				else
					∂R∂Ψ[iZ], ∂R∂Ψ△[iZ], ∂R∂Ψ▽[iZ] = residual.∂RESIDUAL∂Ψ(∂K∂Ψ, discret, hydro, iT, iZ, NiZ, option, param, ΔT, θ, Ψ)
				end # if option.hyPix.∂R∂Ψ_Numerical
			end #for iZ= 1:NiZ

			# FOR TESTING...
				# println("One:=================")
				# println("∂R∂Ψ_Deriv=" , ∂R∂Ψ[1:NiZ],"\n") # No good at cell N
				# println("∂R∂Ψ_Num=" , ∂R∂Ψ2[1:NiZ],"\n")


				# println("Two: =================")
				# println("∂R∂Ψ▽_Num=" , ∂R∂Ψ▽[1:NiZ],"\n")
				# println("∂R∂Ψ▽_Der=" , ∂R∂Ψ▽2[1:NiZ],"\n") # No good

				# println("Tree: =================")
				# println("∂R∂Ψ△_Num=" , ∂R∂Ψ△[1:NiZ],"\n") # Good
				# println("∂R∂Ψ△_Der=" , ∂R∂Ψ△2[1:NiZ],"\n")

		return ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, ΔHpond, θ
		end # function RICHARD
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RESIDUAL_MAX
	#     Averaging the Residuals, depending on method
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RESIDUAL_MAX(discret, iT::Int64, NiZ::Int64, option, Residual, ΔT)
			Residual_Norm = 0.0
			Residual_Max  = 0.0

			# Does not take into consideration the last cell which has a perfect mass balance
			for iZ = 1:NiZ-1
				if option.hyPix.NormMin⍰ == "Norm"
					Residual_Norm += (Residual[iZ] / (ΔT[iT] * discret.ΔZ[iZ])) ^ 2
				else
					Residual_Max = max( Residual_Max, abs(Residual[iZ]) / (ΔT[iT] * discret.ΔZ[iZ]) ) 
				end
			end # for: iZ=NiZ

			if option.hyPix.NormMin⍰ == "Norm"
				Residual_Max = √(Residual_Norm / Float64(NiZ - 1.0))
			end	
		return Residual_Max
		end  # function: RESIDUAL_MAX
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SOLVING_TRIAGONAL_MATRIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT::Int64, iTer::Int64, NiZ::Int64, option, param, Residual, ΔΨmax, θ, Ψ)

			Matrix_Trid = Tridiagonal(∂R∂Ψ△[2:NiZ], ∂R∂Ψ[1:NiZ], ∂R∂Ψ▽[1:NiZ-1])

			Residual = reshape(Residual, NiZ, 1) # Transforming from row to column

			NewtonStep = Matrix_Trid \ -Residual

			for iZ=1:NiZ
				# Iteration k-1
					Ψ₀ = Ψ[iT,iZ]
					θ₀ = θ[iT,iZ]
				
				# Updating Ψ
					if isnan(NewtonStep[iZ])
						@warn isnan(NewtonStep[iZ])
									# println("One:=================")

				# println("iZ = $iZ")

				# println("RESIDUAL =================")
				# println("rESIDUAL=" ,Residual,"\n") # No good at cell N

				
				# println("h: =================")
				# println("Ψ_Deriv=" , Ψ[iT-1, 1:NiZ],"\n") # No good at cell N

				# println("One: =================")
				# println("∂R∂Ψ_Deriv=" , ∂R∂Ψ[1:NiZ],"\n") # No good at cell N

				# println("Two: =================")
				# println("∂R∂Ψ▽_Num=" , ∂R∂Ψ▽[1:NiZ],"\n")

				# println("Tree: =================")
				# println("∂R∂Ψ△_Num=" , ∂R∂Ψ△[1:NiZ],"\n") # Good

						Ψ[iT,iZ] = Ψ₀
					
					else
						Ψ[iT,iZ] += NewtonStep[iZ]

						# Making sure it is within the feasible band 
							Ψ[iT,iZ] = min(max(Ψ[iT,iZ],  param.hyPix.Ψ_MinMin), param.hyPix.Ψ_MaxMax)

						# Correction of θ entering a dry soil 
						if option.hyPix.ZhaWetingDrySoil
							Ψ = ZHA_WETING_DRYSOIL(hydro, iT, iZ, option, θ, θ₀, Ψ, Ψ₀)
						end

						if option.hyPix.DynamicNewtonRaphsonStep
							Ω = DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT, iZ, option, param, ΔΨmax, θ₀, Ψ)
						
							Ψ[iT,iZ] = Ω * Ψ[iT,iZ] + (1.0 - Ω) * Ψ₀
						else
							Ψ[iT,iZ] = param.hyPix.NewtonStep_Mean * Ψ[iT,iZ] + (1.0 - param.hyPix.NewtonStep_Mean) * Ψ₀
						
						end # if option.hyPix.DynamicNewtonRaphsonStep

						if option.hyPix.IterReduceOverShoting
							Ψ = Ψ_REDUCE_OVERSHOOTING(iT, iZ, ΔΨmax, Ψ, Ψ₀)
						end

					end
			end # for iZ=1:NiZ	

		return Ψ
		end  # function: SOLVING_TRIAGONAL_MATRIX
	# ------------------------------------------------------------------

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_REDUCE_OVERSHOOTING
	# 		Making sure that the steps of NR are not too big and within the limits of ΔΨmax
	# 		Does not work
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_REDUCE_OVERSHOOTING(iT, iZ, ΔΨmax, Ψ, Ψold)
			if Ψ[iT,iZ] ≤ Ψold
				Ψ[iT,iZ] = max(Ψ[iT,iZ],  Ψold - ΔΨmax[iZ])
			else
				Ψ[iT,iZ] = min(Ψ[iT,iZ], Ψold + ΔΨmax[iZ])
			end	
		return Ψ
		end  # function: Ψ_REDUCE_OVERSHOOTING
	#---------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NEWTO_NRAPHSON_STEP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT::Int64, iZ::Int64, option, param, ΔΨmax, θ₀, Ψ; Power=1.0)

			θ₁ = Ψ_2_θDual(option.hyPix, Ψ[iT,iZ], iZ, hydro)

			Δθ = abs(θ₁ - θ₀)

			Δθₘₐₓ = timeStep.ΔθMAX(hydro, iT, iZ, option, ΔΨmax, Ψ)
			
		return Ω = param.hyPix.NewtonStep_Max - (param.hyPix.NewtonStep_Max - param.hyPix.NewtonStep_Min) * min(Δθ / Δθₘₐₓ , 1.0) ^ param.hyPix.NewtonStep_Power
		end  # function: NEWTO_NRAPHSON_STEP
	# ------------------------------------------------------------------
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OVERSHOTTING_WET_DRY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	"""Zha, Y., Yang, J., Yin, L., Zhang, Y., Zeng, W., Shi, L., 2017. A modified Picard iteration scheme for overcoming numerical difficulties of simulating infiltration into dry soil. Journal of Hydrology 551, 56–69. https://doi.org/10.1016/j.jhydrol.2017.05.053 """
			function ZHA_WETING_DRYSOIL(hydro, iT, iZ, option, θ, θ₀, Ψ, Ψ₀)
				Ψwet = max( 3.5391 * hydro.σ[iZ]^3 - 20.676 * hydro.σ[iZ]^2 + 24.835 * hydro.σ[iZ] + 15.976, 0.0 )

				Ψdry = exp( 1.6216 * log(hydro.σ[iZ]) + 8.7268 )

				# Determine if there is any oscilation at the wet or dry end of the θ(Ψ) curve
				if  Ψ₀ ≥ Ψdry && Ψ[iT,iZ] ≤ Ψwet
					θ[iT,iZ] = θ₀ + (Ψ[iT,iZ] - Ψ₀) * ∂θ∂Ψ(option.hyPix, Ψ₀, iZ, hydro)

					θ[iT,iZ] = max(min(θ[iT,iZ], hydro.θs[iZ]),  hydro.θr[iZ])

					Ψ[iT,iZ] = θ_2_ΨDual(option.hyPix, θ[iT,iZ] , iZ, hydro)
				end  # Ψ₀ ≥ Ψdry && Ψ[iT,iZ] ≤ Ψwet	
			return Ψ
			end  # function:ZHA_WETING_DRYSOIL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RERUN_HYPIX
	# 		WITH UPDATED Ψ
	#     Rerun if updated ΔT is smaller compared to previously Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RERUN_HYPIX(Count_ReRun::Int64, discret, Flag_NoConverge::Bool, hydro, iT::Int64, NiZ::Int64, option, param, Q, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

			if option.hyPix.Flag_ReRun	&& Count_ReRun ≤ 2
				
				Q[iT,1] = flux.Q!(option, discret, hydro, 1, iT, NiZ, param, ΔHpond, ΔPr, ΔT, θ, Ψ[iT,1], Ψ[iT,1])
				for iZ=1:NiZ
					Q[iT,iZ+1] = flux.Q!(option, discret, hydro, iZ+1, iT, NiZ, param, ΔHpond, ΔPr, ΔT, θ, Ψ[iT, min(iZ+1, NiZ)], Ψ[iT,iZ])
				end

				# for iZ = 1:NiZ
				# 	Q[iT,iZ] = flux.Q!(option, discret, hydro, iZ, iT, NiZ, param, ΔHpond, ΔPr, ΔT, θ, Ψ[iT,iZ], Ψ[iT, max(iZ-1,1)])
				# end # for: iZ= 1:NiZ+1

				ΔT_New, ~ = timeStep.ADAPTIVE_TIMESTEP(discret, hydro, iT, NiZ, option, param, Q, ΔΨmax, ΔSink, θ, Ψ)

				if ΔT[iT] / ΔT_New ≥ param.hyPix.ΔT_Rerun # <>=<>=<>=<>=<>
					ΔT[iT] = ΔT_New
					Flag_ReRun = true
					Count_ReRun += 1

				elseif option.hyPix.TopBoundary⍰ == "Ψ"
					Flag_ReRun = true
					Count_ReRun += 1
				
				else # <>=<>=<>=<>=<>
					Flag_ReRun = false
					Count_ReRun = 0

				end
			else
				Flag_ReRun = false
				Count_ReRun = 0

			end  # if: param.hyPix.ΔT_Rerun

		return Count_ReRun, Flag_ReRun, ΔT
		end  # function: RERUN_HYPIX
	# ------------------------------------------------------------------

end # module: richard
#......................................................................	