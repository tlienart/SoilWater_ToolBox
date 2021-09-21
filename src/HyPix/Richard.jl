# =============================================================
#		MODULE: residual
# =============================================================
module richard
	import ..timeStep, ..flux, ..kunsat, ..ponding, ..residual, ..wrc
	using LinearAlgebra

	export RICHARD_ITERATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD_ITERATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Count_ReRun, discret, Flag_NoConverge::Bool, hydro, iNonConverge::Int64, iT::Int64, Iter_CountTotal::Int64, N_iZ::Int64, param, Q, Residual, Sorptivity::Float64, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, Δθ_Max::Float64, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest, option, optionₘ)

			# INITIALIZING
				for iZ=1:N_iZ
					Ψ[iT,iZ] = Ψbest[iZ]
				end
		
				# Residual_Max_Best it should be improved compared to Ψ[iT,iZ] = Ψ[iT-1,iZ]
				~, ~, ~, ~, Residual, ~, ~ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, hydro, iT, N_iZ, param, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest, option, optionₘ)

				# Averaging the Residuals, depending on method
					Residual_Max_Best = RESIDUAL_MAX(discret, iT, N_iZ, option, Residual, ΔT)
	
			# ITTERATION
			iTer = 0::Int64
			while iTer ≤ param.hyPix.N_Iter	
				iTer += 1
				Iter_CountTotal += 1 # Counting the iterations
		
				∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, ΔHpond, θ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, hydro, iT, N_iZ, param, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest, option, optionₘ)

				Ψ = SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT, iTer, N_iZ, option, param, Residual, ΔΨmax, θ, Ψ)

				# Averaging the Residuals, depending on method
					Residual_Max = RESIDUAL_MAX(discret, iT, N_iZ, option, Residual, ΔT)

				# Determine if itteration made improvement
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
				Count_ReRun, Flag_ReRun, ΔT = RERUN_HYPIX(Count_ReRun, discret, Flag_NoConverge, hydro, iT, N_iZ, option, optionₘ, param, Q, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

		return Count_ReRun, Flag_NoConverge, Flag_ReRun, iNonConverge, Iter_CountTotal, Q, ΔHpond, ΔT, θ, Ψ
		end  # function: RICHARD_SOLVING



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge::Bool, hydro, iT::Int64, N_iZ::Int64, param, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest, option, optionₘ)

			ΔHpond = ponding.PONDING_SORPTIVITY!(discret, hydro, iT, param, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, optionₘ)

			# ∂R∂Ψ2 = fill(0.0, N_iZ)
			# ∂R∂Ψ▽2 = fill(0.0, N_iZ)
			# ∂R∂Ψ△2 =  fill(0.0, N_iZ)

			for iZ=1:N_iZ		
				Q, Residual, θ = residual.RESIDUAL(option, optionₘ, discret, hydro, iT, iZ, N_iZ, param, Q, Residual, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min)

				if option.hyPix.∂R∂Ψ_Numerical				
					∂R∂Ψ[iZ] = residual.∂R∂Ψ_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,N_iZ)], Ψ[iT, iZ], Ψ[iT, min(iZ+1,N_iZ)], Ψ_Max)

					∂R∂Ψ▽[iZ]  = residual.∂R∂Ψ▽_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,N_iZ)], Ψ[iT, iZ], Ψ[iT, min(iZ+1,N_iZ)], Ψ_Max)
		
					∂R∂Ψ△[iZ]  = residual.∂R∂Ψ△_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,N_iZ)], Ψ[iT, iZ], Ψ[iT, min(iZ+1,N_iZ)], Ψ_Max)
				else		
					∂R∂Ψ[iZ], ∂R∂Ψ△[iZ], ∂R∂Ψ▽[iZ] = residual.∂RESIDUAL∂Ψ(∂K∂Ψ, discret, hydro, iT, iZ, N_iZ, option, optionₘ, param, ΔT, θ, Ψ)

				end # if option.hyPix.∂R∂Ψ_Numerical

			end #for iZ= 1:N_iZ

			# println("∂R∂Ψ=" , ∂R∂Ψ[N_iZ]," , ", ∂R∂Ψ2[N_iZ],"\n")
			# println("∂R∂Ψ▽=", ∂R∂Ψ▽[N_iZ] .- ∂R∂Ψ▽2[N_iZ], "\n")
			# println("\n")

			return ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, ΔHpond, θ
		end # function RICHARD

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
				
				# Updating Ψ
					if isnan(NewtonStep[iZ])
						@warn isnan(NewtonStep[iZ])
						Ψ[iT,iZ] = Ψ₀
					
					else
						Ψ[iT,iZ] += NewtonStep[iZ]

						# Making sure it is within the feasible band 
							Ψ[iT,iZ] = min(max(Ψ[iT,iZ], param.hyPix.Ψ_MinMin), param.hyPix.Ψ_MaxMax)

						if option.hyPix.DynamicNewtonRaphsonStep
							θ₀ = θ[iT, iZ]

							θ, Ω = DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT, iZ, option, param, ΔΨmax, θ, θ₀, Ψ)
						
							Ψ[iT,iZ] = Ω * Ψ[iT,iZ] + (1.0 - Ω) * Ψ₀

						else
							Ψ[iT,iZ] = param.hyPix.NewtonStep_Mean * Ψ[iT,iZ] + (1.0 - param.hyPix.NewtonStep_Mean) * Ψ₀
						end # if option.hyPix.DynamicNewtonRaphsonStep

						if option.hyPix.ReduceOvershooting && iTer ≤ 3
							Ψ = Ψ_REDUCE_OVERSHOOTING(iT, iZ, ΔΨmax, Ψ, Ψ₀)
						end
					end
			end # for iZ=1:N_iZ
			
		return Ψ
		end  # function: SOLVING_TRIAGONAL_MATRIX


		
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
		function DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT::Int64, iZ::Int64, option, param, ΔΨmax, θ, θ₀, Ψ)

			θ[iT, iZ] = wrc.Ψ_2_θDual(option.hyPix, Ψ[iT,iZ], iZ, hydro)

			Δθ = abs(θ[iT, iZ] - θ₀)

			Δθₘₐₓ =  timeStep.ΔθMAX(hydro, iT, iZ, option, ΔΨmax, Ψ)
			
			Ω = param.hyPix.NewtonStep_Max - (param.hyPix.NewtonStep_Max - param.hyPix.NewtonStep_Min) * min(Δθ / Δθₘₐₓ , 1.0)
		return θ, Ω
		end  # function: NEWTO_NRAPHSON_STEP
	# ------------------------------------------------------------------
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OVERSHOTTING_WET_DRY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OVERSHOTTING_WET_DRY(hydro, option, iZ, iT)
			
		return
		end  # function: OVERSHOTTING_WET_DRY
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RERUN_HYPIX
	# 		WITH UPDATED Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RERUN_HYPIX(Count_ReRun::Int64, discret, Flag_NoConverge::Bool, hydro, iT::Int64, N_iZ::Int64, option, optionₘ, param, Q, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)
			# Rerun if updated ΔT is smaller compared to previously Ψ

			if option.hyPix.Flag_ReRun	&& Count_ReRun ≤ 3	

				for iZ= 1:N_iZ
					Q[iT,iZ] = flux.Q!(option, optionₘ, discret, hydro, iZ, iT, N_iZ, param, ΔHpond, ΔPr, ΔT, Ψ[iT,iZ], Ψ[iT,max(iZ-1,1)])
				end # for: iZ= 1:N_iZ+1

				ΔT_New, ~ = timeStep.ADAPTIVE_TIMESTEP(discret, hydro, iT, N_iZ, option, option.hyPix, param, Q, ΔΨmax, ΔSink, θ, Ψ)

				# Rerun if the new time step is smaller that the older time step
				# if Flag_NoConverge # <>=<>=<>=<>=<>
				# 	ΔT[iT] = ΔT_New
				# 	Flag_ReRun = true
				# 	Count_ReRun += 1

				if ΔT[iT] / ΔT_New > param.hyPix.ΔT_Rerun # <>=<>=<>=<>=<>
					ΔT[iT] = ΔT_New
					Flag_ReRun = true
					Count_ReRun += 1
				
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

end # module: richard	