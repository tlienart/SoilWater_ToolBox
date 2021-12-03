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
		function RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Count_ReRun::Int64, discret, Flag_NoConverge::Bool, hydro, iNonConverge::Int64, iT::Int64, IterCount::Int64, NiZ::Int64, param, Q, Residual, Sorptivity::Float64, ΔHpond::Vector{Float64}, ΔLnΨmax::Vector{Float64}, ΔPr::Vector{Float64}, ΔSink::Matrix{Float64}, ΔT, θ::Matrix{Float64}, Ψ::Matrix{Float64}, Ψ_Max::Vector{Float64}, Ψbest::Vector{Float64}, option)

			# INITIALIZING
			@inbounds @simd for iZ = 1:NiZ
					Ψ[iT,iZ] = Ψ[iT-1,iZ]
				end # for iZ = 1:NiZ
	
			# ITTERATION
			Residual_Max_Best = Inf
			iTer = 0::Int64
			while iTer ≤ param.hyPix.N_Iter	
				iTer += 1
				IterCount += 1 # Counting the iterations

				if maximum(isnan.(Ψ[iT,1:NiZ]))
					println( iT," , " ,Ψ[iT,1:NiZ])
					error("Ψ=NaN")
				end

				# RESIDUAL MAX BEST: Deriving the Residual max because may be Ψ[iT-1,iZ] is the best solution
				∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, ΔHpond, θ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, hydro, iT, iTer, NiZ, option, param, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max)

				# Computing Residual_Max_Best
				if iTer == 1
					Residual_Max_Best = CONVERGENCECRITERIA(discret, iT, NiZ, option, Residual, ΔT)
				end

				Ψ = SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT, iTer, NiZ, option, param, Residual, ΔLnΨmax, θ, Ψ, Ψ_Max)

				# Averaging the Residuals, depending on method
					Residual_Max = CONVERGENCECRITERIA(discret, iT, NiZ, option, Residual, ΔT)

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

			for iZ=1:NiZ
				θ[iT,iZ] = Ψ_2_θDual(option.hyPix, Ψ[iT,iZ], iZ, hydro)
			end

			# Determine if the simulation is going to rerun with a different time step
				Count_ReRun, Flag_ReRun, ΔT = RERUN_HYPIX(Count_ReRun, discret, Flag_NoConverge, hydro, iT, NiZ, option, param, Q, ΔHpond, ΔLnΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

		return Count_ReRun, Flag_NoConverge, Flag_ReRun, iNonConverge, IterCount, Q, ΔHpond, ΔT, θ, Ψ
		end  # function: RICHARD_SOLVING
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge::Bool, hydro, iT::Int64, iTer::Int64, NiZ::Int64, option, param, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max)

			if option.hyPix.Ponding
				ΔHpond = ponding.PONDING_SORPTIVITY!(discret, hydro, iT, option, param, Q, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ)
			else
				ΔHpond[iT] = 0.0
			end

	#----------------------------------------------------------------
			# ∂R∂Ψ2 = fill(0.0, NiZ)
			# ∂R∂Ψ▽2 = fill(0.0, NiZ)
			# ∂R∂Ψ△2 =  fill(0.0, NiZ)

			for iZ=1:NiZ
				Q, Residual, θ = residual.RESIDUAL(discret, hydro, iT, iZ, NiZ, option, param, Q, Residual, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ)

				if option.hyPix.∂R∂Ψ_Numerical
					∂R∂Ψ[iZ] = residual.∂R∂Ψ_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, option, param, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,NiZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,NiZ)], Ψ_Max)

					∂R∂Ψ▽[iZ]  = residual.∂R∂Ψ▽_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, option, param, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,NiZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,NiZ)], Ψ_Max)
		
					∂R∂Ψ△[iZ]  = residual.∂R∂Ψ△_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, option, param, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,NiZ)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,NiZ)], Ψ_Max)
				else
					∂R∂Ψ[iZ], ∂R∂Ψ△[iZ], ∂R∂Ψ▽[iZ] = residual.∂RESIDUAL∂Ψ(∂K∂Ψ, discret, hydro, iT, iZ, NiZ, option, param, ΔT, θ, Ψ)
				end # if option.hyPix.∂R∂Ψ_Numerical
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

		return ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, ΔHpond, θ
		end # function RICHARD
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERGENCECRITERIA
	#     Averaging the Residuals, depending on method
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CONVERGENCECRITERIA(discret, iT::Int64, NiZ::Int64, option, Residual, ΔT)
			Residual_Norm = 0.0
			Residual_Max  = 0.0

			# Does not take into consideration the last cell which has a perfect mass balance
			for iZ = 1:NiZ-1
				if option.hyPix.NormMin⍰ == "Norm"
					Residual_Norm += (Residual[iZ] / (ΔT[iT] * discret.ΔZ[iZ])) ^ 2.0
				else
					Residual_Max = max( Residual_Max, abs(Residual[iZ]) / (ΔT[iT] * discret.ΔZ[iZ]) ) 
				end
			end # for: iZ=NiZ

			if option.hyPix.NormMin⍰ == "Norm"
				return Residual_Max = √(Residual_Norm / Float64(NiZ - 1.0))
			end
					
		return Residual_Max
		end  # function: CONVERGENCECRITERIA
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SOLVING_TRIAGONAL_MATRIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT::Int64, iTer::Int64, NiZ::Int64, option, param, Residual, ΔLnΨmax, θ, Ψ, Ψ_Max)

			# if maximum(isnan.(Ψ[iT, 1:NiZ]))
			# 	println("Ψ=$(Ψ[iT, 1:NiZ]) \n")
			# 	println("0 Pressure")
			# end

			# if maximum(isnan.(∂R∂Ψ△[2:NiZ]))
			# 	println("Ψ=$(Ψ[iT, 2:NiZ]) \n")
			# 	println(∂R∂Ψ△[2:NiZ])
			# 	println("1 ∂R∂Ψ△")
			# end

			# if maximum(isnan.(∂R∂Ψ[1:NiZ]))
			# 	println("Ψ=$(Ψ[iT, 1:NiZ]) \n")
			# 	println(∂R∂Ψ[1:NiZ])
			# 	println("2 ∂R∂Ψ")
			# end

			# if maximum(isnan.(∂R∂Ψ▽[1:NiZ-1]))
			# 	println("Ψ=$(Ψ[iT, 1:NiZ]) \n")
			# 	println(∂R∂Ψ▽[2:NiZ])
			# 	println("3 ∂R∂Ψ▽")
			# end

			Matrix_Trid = Tridiagonal(∂R∂Ψ△[2:NiZ], ∂R∂Ψ[1:NiZ], ∂R∂Ψ▽[1:NiZ-1])

			# if maximum(isnan.(Matrix_Trid))
			# 	println("Ψ=$(Ψ[iT, 1:NiZ]) \n")
			# 	println(Matrix_Trid)
			# 	println("4 Tridiagonal")
			# end

			Residual = reshape(Residual, NiZ, 1) # Transforming from row to column

			# if maximum(isnan.(Residual))
			# 	println("Ψ=$(Ψ[iT, 1:NiZ]) \n")
			# 	println(Residual)
			# 	println("5 Residual")
			# end

			NewtonStep = Matrix_Trid \ -Residual

			# # More robust method for computing the Newton step if the previous fails
			# 	if maximum(isnan.(NewtonStep))
			# 		NewtonStep = lu(Matrix_Trid) \ -Residual
			# 	end

			# if maximum(isnan.(NewtonStep))
			# 	println("Ψ=$(Ψ[iT, 1:NiZ]) \n")
			# 	println(NewtonStep)
			# 	println("6 NewtonStep")
			# end

			# if isnan(Ψ[iT,iZ])
			# 	println( iT," , " ,Ψ[iT,iZ])
			# 	error("Richard Ψ=NaN")
			# end

			for iZ=1:NiZ
				# Iteration k-1
					Ψ₀ = Ψ[iT,iZ]
					θ₀ = θ[iT,iZ]
				
				# Updating Ψ
				if isnan(NewtonStep[iZ])

					# println("iZ = $iZ")
					# println("NewtonStep = $(NewtonStep[iZ]), \n")
					# println("h: =================")
					# println("Ψ=" , Ψ[iT, 1:NiZ],"\n")
					# println("RESIDUAL =================")
					# println("RESIDUAL=" ,Residual,"\n")				
					# println("One: =================")
					# println("∂R∂Ψ_Deriv=" , ∂R∂Ψ[1:NiZ],"\n") # No good at cell N
					# println("Two: =================")
					# println("∂R∂Ψ▽_Num=" , ∂R∂Ψ▽[1:NiZ],"\n")
					# println("Tree: =================")
					# println("∂R∂Ψ△_Num=" , ∂R∂Ψ△[1:NiZ],"\n") # Good
					
					Ψ[iT,iZ] = Ψ₀ + eps(100.0)
					
					
				else
					# Newtyon step
						Ψ[iT,iZ] += NewtonStep[iZ]

					# Correction of θ entering a dry soil 
					if option.hyPix.ZhaWetingDrySoil
						Ψ = ZHA_WETING_DRYSOIL(hydro, iT, iZ, option, θ, θ₀, Ψ, Ψ₀)
					end

					# Making sure it is within the feasible band 
						Ψ[iT,iZ] = min(max(Ψ[iT,iZ], param.hyPix.Ψ_MinMin), Ψ_Max[iZ])

					if option.hyPix.DynamicNewtonRaphsonStep
						Ω = DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT, iZ, option, param, ΔLnΨmax, θ₀, Ψ)

						# if isnan(Ω)
						# 	println(Ω)
						# 	error("DYNAMIC_NEWTON_RAPHSON_STEP")
						# end
					
						Ψ[iT,iZ] = Ω * Ψ[iT,iZ] + (1.0 - Ω) * Ψ₀
						
					else
						Ψ[iT,iZ] = 0.5 * (Ψ[iT,iZ] +  Ψ₀)

					end # if option.hyPix.DynamicNewtonRaphsonStep

					if option.hyPix.IterReduceOverShoting
						Ψ = Ψ_REDUCE_OVERSHOOTING(iT, iZ, ΔLnΨmax, Ψ, Ψ₀)
					end

				end
			end # for iZ=1:NiZ	

		return Ψ
		end  # function: SOLVING_TRIAGONAL_MATRIX
	# ------------------------------------------------------------------

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_REDUCE_OVERSHOOTING
	# 		Making sure that the steps of NR are not too big and within the limits of ΔLnΨmax
	# 		Does not work
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_REDUCE_OVERSHOOTING(iT, iZ, ΔLnΨmax, Ψ, Ψ₀)
			if Ψ[iT,iZ] ≤ Ψ₀
				Ψ[iT,iZ] = max(Ψ[iT,iZ], expm1(log1p(Ψ₀) - ΔLnΨmax[iZ]))
			else
				Ψ[iT,iZ] = min(Ψ[iT,iZ], expm1(log1p(Ψ₀) + ΔLnΨmax[iZ]))
			end	
		return Ψ
		end  # function: Ψ_REDUCE_OVERSHOOTING
	#---------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NEWTO_NRAPHSON_STEP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT::Int64, iZ::Int64, option, param, ΔLnΨmax, θ₀, Ψ)

			θ₁ = Ψ_2_θDual(option.hyPix, Ψ[iT,iZ], iZ, hydro)

			Δθ = abs(θ₁ - θ₀)

			Δθₘₐₓ = timeStep.ΔθMAX(hydro, iT, iZ, option, ΔLnΨmax, Ψ) 
			
		return param.hyPix.NewtonStep_Max - (param.hyPix.NewtonStep_Max - param.hyPix.NewtonStep_Min) * min(Δθ / Δθₘₐₓ, 1.0) ^ param.hyPix.NewtonStep_Power
		end  # function: NEWTO_NRAPHSON_STEP
	# ------------------------------------------------------------------
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OVERSHOTTING_WET_DRY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	"""Zha, Y., Yang, J., Yin, L., Zhang, Y., Zeng, W., Shi, L., 2017. A modified Picard iteration scheme for overcoming numerical difficulties of simulating infiltration into dry soil. Journal of Hydrology 551, 56–69. https://doi.org/10.1016/j.jhydrol.2017.05.053 """
			function ZHA_WETING_DRYSOIL(hydro, iT, iZ, option, θ, θ₀, Ψ, Ψ₀)
				# Ψwet = max( 3.5391 * hydro.σ[iZ]^3 - 20.676 * hydro.σ[iZ]^2 + 24.835 * hydro.σ[iZ] + 15.976, 0.0 )

				Ψwet = max(-2.3116 * hydro.σ[iZ] ^ 2.0 - 2.9372 * hydro.σ[iZ] + 27.83, 0.0)

				Ψdry = exp(1.6216 * log(hydro.σ[iZ]) + 8.7268)

				# Determine if there is any oscilation at the wet or dry end of the θ(Ψ) curve
				if (Ψ[iT,iZ] ≤ Ψwet && Ψ₀ ≥ Ψdry)
					θ[iT,iZ] = θ₀ + (Ψ[iT,iZ] - Ψ₀) * ∂θ∂Ψ(option.hyPix, Ψ₀, iZ, hydro)

					θ[iT,iZ] = max(min(θ[iT,iZ], hydro.θs[iZ]), hydro.θr[iZ])

					Ψ[iT,iZ] = θ_2_ΨDual(option.hyPix, θ[iT,iZ] , iZ, hydro)
				end  # Ψ[iT,iZ] ≤ Ψwet && Ψ₀ ≥ Ψdry

				# if isnan(Ψ[iT,iZ])
				# 	error("ZHA_WETING_DRYSOIL Ψ= NaN, iT= $iT, iZ= $iZ")
				# end

			return Ψ
			end  # function:ZHA_WETING_DRYSOIL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RERUN_HYPIX
	# 		WITH UPDATED Ψ
	#     Rerun if updated ΔT is smaller compared to previously Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RERUN_HYPIX(Count_ReRun::Int64, discret, Flag_NoConverge::Bool, hydro, iT::Int64, NiZ::Int64, option, param, Q, ΔHpond, ΔLnΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

			if option.hyPix.Flag_ReRun	&& Count_ReRun ≤ 3		
				Q[iT,1] = flux.Q!(option, discret, hydro, 1, iT, NiZ, param, Q, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT,1], Ψ[iT,1])
				for iZ=1:NiZ
					Q[iT,iZ+1] = flux.Q!(option, discret, hydro, iZ+1, iT, NiZ, param, Q, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, min(iZ+1, NiZ)], Ψ[iT,iZ])
				end

				ΔT_New, ~ = timeStep.ADAPTIVE_TIMESTEP(discret, hydro, iT, NiZ, option, param, Q, ΔLnΨmax, ΔSink, θ, Ψ)

				if ΔT[iT] / ΔT_New ≥ param.hyPix.ΔT_Rerun # <>=<>=<>=<>=<>
					ΔT[iT] = ΔT_New
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