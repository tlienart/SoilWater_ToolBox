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
		function RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Count_ReRun::Int64, discret, Flag_NoConverge::Bool, hydro, iNonConverge::Int64, iT::Int64, IterCount::Int64, N_iZ::Int64, param, Q, Residual, Sorptivity::Float64, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, Δθ_Max::Float64, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest, option)

			# INITIALIZING
				for iZ=1:N_iZ
					if option.hyPix.TopBoundary⍰ == "Ψ"
						Ψbest[1] =  param.hyPix.Ψ_Top
					end
					Ψ[iT,iZ] = Ψbest[iZ]
				end
		
				# Residual_Max_Best it should be improved compared to Ψ[iT,iZ] = Ψ[iT-1,iZ]
				~, ~, ~, ~, Residual, ~, ~ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, hydro, iT, N_iZ, param, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest, option)

				# Averaging the Residuals, depending on method
					Residual_Max_Best = RESIDUAL_MAX(discret, iT, N_iZ, option, Residual, ΔT)
	
			# ITTERATION
			iTer = 0::Int64
			while iTer ≤ param.hyPix.N_Iter	
				iTer += 1
				IterCount += 1 # Counting the iterations
		
				∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, ΔHpond, θ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, hydro, iT, N_iZ, param, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest, option)

				Ψ = SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT, iTer, N_iZ, option, param, Residual, ΔΨmax, θ, Ψ)

				# Averaging the Residuals, depending on method
					Residual_Max = RESIDUAL_MAX(discret, iT, N_iZ, option, Residual, ΔT)

				# Determine if iteration made improvement
					if Residual_Max < Residual_Max_Best	
						for iZ=1:N_iZ
							Ψbest[iZ] = Ψ[iT,iZ]
						end
						Residual_Max_Best = Residual_Max
					end # Residual_Max < Residual_Max_Best 	

				# Did we achieve the goals
				if Residual_Max < param.hyPix.WaterBalanceResidual_Max
					break # Move out the loop
				end  # if: Residual
			end # while: iTer ======================

			# Making sure we get the best if convergence fails
			if iTer ≥ param.hyPix.N_Iter + 1
				iNonConverge += 1

				# if non converge compute Q(Ψbest)
				# if option.hyPix.NoConverge_Ψbest
				# 	Flag_NoConverge = true
				# end

				# Put the best values
				for iZ=1:N_iZ
					Ψ[iT,iZ] = Ψbest[iZ]
				end
			else
				Flag_NoConverge = false
			end #  iTer == param.hyPix.N_Iter

			if option.hyPix.Qbottom_Correction
				Q[iT,N_iZ+1] = max(( - ΔSink[iT,N_iZ] - discret.ΔZ[N_iZ] * ((θ[iT,N_iZ] - θ[iT-1,N_iZ]) - hydro.So[N_iZ] * (Ψ[iT,N_iZ]- Ψ[iT-1,N_iZ]) * (θ[iT,N_iZ] / hydro.θs[N_iZ]))) / ΔT[iT] + Q[iT,N_iZ], 0.0)
			end

			# Determine if the simulation is going to rerun with a different time step
				Count_ReRun, Flag_ReRun, ΔT = RERUN_HYPIX(Count_ReRun, discret, Flag_NoConverge, hydro, iT, N_iZ, option, param, Q, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

		return Count_ReRun, Flag_NoConverge, Flag_ReRun, iNonConverge, IterCount, Q, ΔHpond, ΔT, θ, Ψ
		end  # function: RICHARD_SOLVING
	#-----------------------------------------------------------------



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge::Bool, hydro, iT::Int64, N_iZ::Int64, param, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest, option)

			ΔHpond = ponding.PONDING_SORPTIVITY!(discret, hydro, iT, param, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, option)

	#----------------------------------------------------------------

			# ∂R∂Ψ2 = fill(0.0, N_iZ)
			# ∂R∂Ψ▽2 = fill(0.0, N_iZ)
			# ∂R∂Ψ△2 =  fill(0.0, N_iZ)

			for iZ=1:N_iZ		
				Q, Residual, θ = residual.RESIDUAL(option, discret, hydro, iT, iZ, N_iZ, param, Q, Residual, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min)

				if option.hyPix.∂R∂Ψ_Numerical				
					∂R∂Ψ[iZ] = residual.∂R∂Ψ_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, option, param, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,N_iZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,N_iZ)], Ψ_Max)

					∂R∂Ψ▽[iZ]  = residual.∂R∂Ψ▽_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, option, param, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,N_iZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,N_iZ)], Ψ_Max)
		
					∂R∂Ψ△[iZ]  = residual.∂R∂Ψ△_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, option, param, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,N_iZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,N_iZ)], Ψ_Max)
				else		
					∂R∂Ψ[iZ], ∂R∂Ψ△[iZ], ∂R∂Ψ▽[iZ] = residual.∂RESIDUAL∂Ψ(∂K∂Ψ, discret, hydro, iT, iZ, N_iZ, option, param, ΔT, θ, Ψ)

				end # if option.hyPix.∂R∂Ψ_Numerical

			end #for iZ= 1:N_iZ

			# println("∂R∂Ψ=" , ∂R∂Ψ[N_iZ]," , ", ∂R∂Ψ2[N_iZ],"\n")
			# println("∂R∂Ψ▽=", ∂R∂Ψ▽[N_iZ] .- ∂R∂Ψ▽2[N_iZ], "\n")
			# println("\n")

		return ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, ΔHpond, θ
		end # function RICHARD
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RESIDUAL_MAX
	#     Averaging the Residuals, depending on method
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RESIDUAL_MAX(discret, iT::Int64, N_iZ::Int64, option, Residual, ΔT)
				Residual_Norm = 0.0
				Residual_Max  = 0.0

				# Does not take into consideration the last cell which has a perfect mass balance
				for iZ=1:N_iZ
					if option.hyPix.NormMin⍰ == "Norm"
						Residual_Norm += (Residual[iZ] / (ΔT[iT] * discret.ΔZ[iZ])) ^ 2
					else
						Residual_Max = max( Residual_Max, abs(Residual[iZ]) / (ΔT[iT] * discret.ΔZ[iZ]) ) 
					end
				end # for: iZ=N_iZ

				if option.hyPix.NormMin⍰ == "Norm"
					Residual_Max = √(Residual_Norm / N_iZ)
				end
			
		return Residual_Max
		end  # function: RESIDUAL_MAX
	#-----------------------------------------------------------------



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SOLVING_TRIAGONAL_MATRIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT::Int64, iTer::Int64, N_iZ::Int64, option, param, Residual, ΔΨmax, θ, Ψ)

			Matrix_Trid = Tridiagonal(∂R∂Ψ△[2:N_iZ], ∂R∂Ψ[1:N_iZ], ∂R∂Ψ▽[1:N_iZ-1])

			Residual = reshape(Residual, N_iZ, 1) # Transforming from row to column

			NewtonStep = Matrix_Trid \ -Residual

			for iZ=1:N_iZ
				# Iteration k-1
					Ψ₀ = Ψ[iT,iZ]
					θ₀ = θ[iT,iZ]
				
				# Updating Ψ
					if isnan(NewtonStep[iZ])
						@warn isnan(NewtonStep[iZ])
						Ψ[iT,iZ] = Ψ₀
					
					else
						Ψ[iT,iZ] += NewtonStep[iZ]

						if option.hyPix.IterReduceOverShoting && iTer ≤ 3
							Ψ = Ψ_REDUCE_OVERSHOOTING(iT, iZ, ΔΨmax, Ψ, Ψ₀)
						end

						# Correction of θ entering a dry soil 
						if option.hyPix.ZhaWetingDrySoil
							Ψ = ZHA_WETING_DRYSOIL(hydro, iT, iZ, option, θ, θ₀, Ψ, Ψ₀)
						end

						# Making sure it is within the feasible band 
							Ψ[iT,iZ] = min(max(Ψ[iT,iZ], eps(10.0)), param.hyPix.Ψ_MaxMax)

						if option.hyPix.DynamicNewtonRaphsonStep
							Ω = DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT, iZ, option, param, ΔΨmax, θ₀, Ψ)
						
							Ψ[iT,iZ] = Ω * Ψ[iT,iZ] + (1.0 - Ω) * Ψ₀

						else
							Ψ[iT,iZ] = param.hyPix.NewtonStep_Mean * Ψ[iT,iZ] + (1.0 - param.hyPix.NewtonStep_Mean) * Ψ₀
						end # if option.hyPix.DynamicNewtonRaphsonStep

					end
			end # for iZ=1:N_iZ	

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
		function DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT::Int64, iZ::Int64, option, param, ΔΨmax, θ₀, Ψ)

			θ₁ = Ψ_2_θDual(option.hyPix, Ψ[iT,iZ], iZ, hydro)

			Δθ = abs(θ₁ - θ₀)

			Δθₘₐₓ = timeStep.ΔθMAX(hydro, iT, iZ, option, ΔΨmax, Ψ)
			
		return Ω = param.hyPix.NewtonStep_Max - (param.hyPix.NewtonStep_Max - param.hyPix.NewtonStep_Min) * min(Δθ / Δθₘₐₓ , 1.0)
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
		function RERUN_HYPIX(Count_ReRun::Int64, discret, Flag_NoConverge::Bool, hydro, iT::Int64, N_iZ::Int64, option, param, Q, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

			if option.hyPix.Flag_ReRun	&& Count_ReRun ≤ 2	
				for iZ = 1:N_iZ
					Q[iT,iZ] = flux.Q!(option, discret, hydro, iZ, iT, N_iZ, param, ΔHpond, ΔPr, ΔT, Ψ[iT,iZ], Ψ[iT, max(iZ-1,1)])
				end # for: iZ= 1:N_iZ+1

				ΔT_New, ~ = timeStep.ADAPTIVE_TIMESTEP(discret, hydro, iT, N_iZ, option, param, Q, ΔΨmax, ΔSink, θ, Ψ)

				if ΔT[iT] / ΔT_New ≥ param.hyPix.ΔT_Rerun # <>=<>=<>=<>=<>
					ΔT[iT] = ΔT_New
					Flag_ReRun = true
					Count_ReRun += 1

				elseif option.hyPix.TopBoundary⍰ == "Ψ"
					Flag_ReRun = true
					Count_ReRun += 2
				
				else # <>=<>=<>=<>=<>
					Flag_ReRun = false
					Flag_NoConverge = false
					Count_ReRun = 0
				end
			else
				Flag_ReRun = false
				Flag_NoConverge = false
				Count_ReRun = 0
			end  # if: param.hyPix.ΔT_Rerun

		return Count_ReRun, Flag_ReRun, ΔT
		end  # function: RERUN_HYPIX
	# ------------------------------------------------------------------

end # module: richard
#......................................................................	