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
		function RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, iCount_ReRun::Int64, discret, Flag_NoConverge::Bool, hydro, iNonConverge::Int64, iT::Int64, IterCount::Int64, NiZ::Int64, paramHypix, Q, Residual, Sorptivity::Float64, Hpond::Vector{Float64}, ΔLnΨmax::Vector{Float64}, ΔPr::Vector{Float64}, ΔRunoff::Vector{Float64}, ΔSink::Matrix{Float64}, ΔT, θ::Matrix{Float64}, Ψ::Matrix{Float64},Ψ_Min::Vector{Float64}, Ψ_Max::Vector{Float64}, Ψbest::Vector{Float64}, optionHypix)
						
			# INITIALIZING
			@inbounds @simd for iZ = 1:NiZ
					Ψ[iT,iZ] = Ψ[iT-1,iZ]
				end # for iZ = 1:NiZ
	
			# ITTERATION
			Residual_Max_Best = Inf
			iTer = 0::Int64
			while iTer ≤ paramHypix.N_Iter - 1	
            iTer      += 1
            IterCount += 1 # Counting the iterations

				# RESIDUAL MAX BEST: Deriving the Residual max because may be Ψ[iT-1,iZ] is the best solution
				∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, Hpond, θ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, hydro, iT, NiZ, optionHypix, paramHypix, Q, Residual, Sorptivity, Hpond, ΔPr, ΔRunoff, ΔSink, ΔT, θ, Ψ, Ψ_Max)

				# Computing Residual_Max_Best at the beginning before iteration
				if iTer == 1
					Residual_Max_Best = CONVERGENCECRITERIA(discret, iT, NiZ, Residual, ΔT)
				end

				# The minimum Ψ depends if we are close to saturation
					Ψ_Min = ΨMIN(iT, NiZ, paramHypix, Ψ, Ψ_Min)

					Ψ = SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT, NiZ, optionHypix, paramHypix, Residual, ΔLnΨmax, θ, Ψ, Ψ_Min, Ψ_Max)

				# Averaging the Residuals, depending on method
					Residual_Max = CONVERGENCECRITERIA(discret, iT, NiZ,  Residual, ΔT)

				# Determine if iteration made improvement
					if Residual_Max < Residual_Max_Best	
						@inbounds @simd for iZ=1:NiZ
							Ψbest[iZ] = Ψ[iT,iZ]
						end
						Residual_Max_Best = copy(Residual_Max)
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

			# UPDATE Θ
			for iZ=1:NiZ
				θ[iT,iZ] = Ψ_2_θDual(optionHypix, Ψ[iT,iZ], iZ, hydro)
			end

			# Determine if the simulation is going to rerun with a different time step
			Flag_ReRun, iCount_ReRun, iNonConverge, ΔT = RERUN_HYPIX(discret, Flag_NoConverge, Hpond, hydro, iCount_ReRun, iNonConverge, iT, NiZ, optionHypix, paramHypix, Q, ΔLnΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

		return Flag_NoConverge, Flag_ReRun, Hpond, iCount_ReRun, iNonConverge, iTer, IterCount, Q, ΔRunoff, ΔT, θ, Ψ
		end  # function: RICHARD_SOLVING
	#----------------------------------------------------------]-------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge::Bool, hydro, iT::Int64, NiZ::Int64, optionHypix, paramHypix, Q, Residual, Sorptivity, Hpond, ΔPr, ΔRunoff, ΔSink, ΔT, θ, Ψ, Ψ_Max)

			if optionHypix.Ponding && (ΔPr[iT] > eps(100.) || Hpond[iT-1] > eps(100.))
				Hpond, ΔRunoff = ponding.PONDING_RUNOFF_SORPTIVITY(discret, Hpond, hydro, iT, optionHypix, paramHypix, Sorptivity, ΔPr, ΔRunoff, ΔSink, ΔT, θ, Ψ)
			else
            Hpond[iT]   = 0.0::Float64
            ΔRunoff[iT] = 0.0::Float64
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
	#		FUNCTION : ΨMIN
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ΨMIN(iT::Int64, NiZ::Int64, paramHypix, Ψ::Matrix{Float64}, Ψ_Min::Vector{Float64})
			@inbounds @simd for iZ=1:NiZ
				if Ψ[iT-1,iZ] < paramHypix.opt.ΨmacMat / 2.0 # mm
					Ψ_Min[iZ] = paramHypix.Ψ_MinMin
				else
					Ψ_Min[iZ] = 0.0
				end
			end
		return Ψ_Min
		end  # function: ΨMIN
	# ------------------------------------------------------------------

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERGENCECRITERIA
	#     Averaging the Residuals, depending on method
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CONVERGENCECRITERIA(discret, iT::Int64, NiZ::Int64, Residual, ΔT)
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
		function SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT::Int64, NiZ::Int64, optionHypix, paramHypix, Residual, ΔLnΨmax, θ, Ψ, Ψ_Min::Vector{Float64}, Ψ_Max::Vector{Float64})

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
				
					# Correction of θ entering a dry soil 
						Ψ = ZHA_WETING_DRYSOIL(hydro, iT, iZ, optionHypix, θ, θ₀, Ψ, Ψ₀)

					# Assuring that the limits of Ψ are physical
						Ψ[iT,iZ] = min(max(Ψ[iT,iZ], Ψ_Min[iZ]), Ψ_Max[iZ])

					# Smootening the steps
						Ω = DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT, iZ, optionHypix, ΔLnΨmax, θ₀, Ψ)		
						Ψ[iT,iZ] = Ω * Ψ[iT,iZ] + (1.0 - Ω) * Ψ₀

				# No comvergence
				else
					# @warn error("===== Difficulties in inverting Tridiagonal =====")
					Ψ[iT,iZ] = Ψ₀ + eps(100.0)
					println(" ================   STRUGGLING ====================")

				end
			end # for iZ=1:NiZ	
		return Ψ
		end  # function: SOLVING_TRIAGONAL_MATRIX
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NEWTO_NRAPHSON_STEP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT::Int64, iZ::Int64, optionHypix, ΔLnΨmax, θ₀, Ψ)

			θ₁ = Ψ_2_θDual(optionHypix, Ψ[iT,iZ], iZ, hydro)

			Δθ = abs(θ₁ - θ₀)

			Δθₘₐₓ = timeStep.ΔθMAX(hydro, iT, iZ, optionHypix, ΔLnΨmax, Ψ) 
		return 1.0 - 0.8 * min(Δθ / Δθₘₐₓ, 1.0) ^ 2.0
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
		function RERUN_HYPIX(discret, Flag_NoConverge::Bool, Hpond::Vector{Float64}, hydro, iCount_ReRun::Int64, iNonConverge::Int64, iT::Int64, NiZ::Int64, optionHypix, paramHypix, Q::Matrix{Float64}, ΔLnΨmax::Vector{Float64}, ΔPr::Vector{Float64}, ΔSink::Matrix{Float64}, ΔT::Vector{Float64}, θ::Matrix{Float64}, Ψ::Matrix{Float64})

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

			if iCount_ReRun ≤ 2	
				ΔTₒ = COMPUTE_ΔT(discret, hydro, iT, NiZ, optionHypix, paramHypix, Q, Hpond, ΔLnΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

				if ΔTₒ < paramHypix.ΔT_Min + paramHypix.ΔT_MaxChange * max(ΔT[iT] - paramHypix.ΔT_Min, 0.0)
               Flag_ReRun     = true
               ΔT[iT]         = ΔTₒ
               iCount_ReRun  += 1
				
				elseif Flag_NoConverge
               Flag_ReRun     = true
               ΔT[iT]         = paramHypix.ΔT_Min + 0.8 * max(ΔT[iT] - paramHypix.ΔT_Min, 0.0)
               iCount_ReRun   += 1

				else # <>=<>=<>=<>=<>
               Flag_ReRun   = false
               iCount_ReRun = 1
				end
			else
				Flag_ReRun = false
				iCount_ReRun = 1

				if Flag_NoConverge
					iNonConverge += 1
				end
			end  # if: paramHypix.ΔT_MaxChange

		return Flag_ReRun, iCount_ReRun, iNonConverge, ΔT
		end  # function: RERUN_HYPIX
	# ------------------------------------------------------------------

end # module: richard
#......................................................................