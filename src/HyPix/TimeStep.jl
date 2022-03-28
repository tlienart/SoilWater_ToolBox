module timeStep
	import ..wrc
   export TIMESTEP, ADAPTIVE_TIMESTEP, ΔΨMAX!

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION :  TIMESTEP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIMESTEP(∑T::Vector{Float64}, discret, Flag_ReRun::Bool, hydro, iT::Int64, iTer::Int64, N_∑T_Climate::Float64, NiZ::Int64, optionHypix, paramHypix, Q::Matrix{Float64}, ΔLnΨmax::Vector{Float64}, ΔSink::Matrix{Float64}, ΔT::Vector{Float64}, θ::Matrix{Float64}, Ψ::Matrix{Float64})

			Δθ_Max = paramHypix.Δθ_Max

			# The iT is of the previous simulation
			if !Flag_ReRun # <>=<>=<>=<>=<>	
				ΔT₂, Δθ_Max = ADAPTIVE_TIMESTEP(discret, hydro, iT, NiZ, optionHypix, paramHypix, Q, ΔLnΨmax, ΔSink, θ, Ψ)
				iT += 1 # Going to the next simulation
				ΔT[iT] = ΔT₂
			end

			# Check if we are at the last time step
			if N_∑T_Climate - (∑T[iT-1] + ΔT[iT]) <= 0.00001
				if N_∑T_Climate - ∑T[iT-1] < 0.00001
					ΔT[iT] = eps()
					FlagContinueLoop = false
				else # New time step
					ΔT[iT] = N_∑T_Climate - ∑T[iT-1]
					∑T[iT] = ∑T[iT-1] + ΔT[iT]
					FlagContinueLoop = true
				end
			else # Not at the last time step: N_∑T_Climate - (∑T[iT] + ΔT) > 0.0
				∑T[iT] = ∑T[iT-1] + ΔT[iT]
				FlagContinueLoop = true
			end #  N_∑T_Climate - (∑T[iT] + ΔT) < 0.0

		return ∑T, FlagContinueLoop, iT, ΔT, Δθ_Max
		end # TIMESTEP()
      

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΔΨMAX!
	# 		Computing ΔΨMAX required by ADAPTIVE_TIMESTEP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ΔΨMAX!(hydro, NiZ::Int64, optionHypix, paramHypix, ΔLnΨmax::Vector{Float64})
			for iZ=1:NiZ
				θ½ = (hydro.θsMacMat[iZ] + hydro.θr[iZ]) * 0.5
				
				θ△ = min(θ½ + paramHypix.Δθ_Max * 0.5, hydro.θs[iZ])

				θ▽ = max(θ½ - paramHypix.Δθ_Max * 0.5, hydro.θr[iZ])

				ΔLnΨmax[iZ] = (log1p(wrc.θ_2_ΨDual(optionHypix, θ▽, iZ, hydro)) - log1p(wrc.θ_2_ΨDual(optionHypix, θ△, iZ, hydro))) * 0.5	
			end # for iZ=1:NiZ	
		return ΔLnΨmax
		end  # function: ΔΨMAX!
	#--------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΔθMAX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ΔθMAX(hydro, iT::Int64, iZ::Int64, optionHypix, ΔLnΨmax::Vector{Float64}, Ψ::Matrix{Float64})
			Ψ₀ = max(Ψ[iT,iZ], 0.0::Float64)

			if log1p(Ψ₀) > ΔLnΨmax[iZ]
				Ψ▽ = expm1(log1p(Ψ₀) - ΔLnΨmax[iZ])		
			else
				Ψ▽ = 0.0::Float64
			end	

			Ψ△  = expm1(log1p(Ψ₀) + ΔLnΨmax[iZ])
			
			θ△ = wrc.Ψ_2_θDual(optionHypix, Ψ▽, iZ, hydro) + eps(100.0)
		
			θ▽ = wrc.Ψ_2_θDual(optionHypix, Ψ△, iZ, hydro)
		return θ△ - θ▽
		end  # function:  ΔθMAX
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ADAPTIVE_TIMESTEP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ADAPTIVE_TIMESTEP(discret, hydro, iT::Int64, NiZ::Int64, optionHypix, paramHypix, Q, ΔLnΨmax, ΔSink, θ, Ψ)
			
			# Initializing
				Δθ₂_Max = paramHypix.Δθ_Max
				ΔT_New_Norm = 0.0::Float64

			# Computing smallest Δθ_Max
				Ngood = 0::Int64
				for iZ = 1:NiZ-1
					if abs(Ψ[iT,iZ] - Ψ[iT,iZ+1]) ≥ 1.0 # mm
						if optionHypix.AdaptiveTimeStep⍰ == "ΔΨ" # <>=<>=<>=<>=<>
							Δθ₂_Max = ΔθMAX(hydro, iT, iZ, optionHypix, ΔLnΨmax, Ψ)	
						end # optionHypix.AdaptiveTimeStep⍰ ==:ΔΨ
					
						ΔT₂_New = (discret.ΔZ[iZ] * Δθ₂_Max + ΔSink[iT,iZ]) / (abs(Q[iT,iZ] - Q[iT,iZ+1]))

						ΔT₂_New = min(max(paramHypix.ΔT_Min, ΔT₂_New), paramHypix.ΔT_Max)
		
						ΔT_New_Norm += ΔT₂_New ^ 2.0

						Ngood += 1
					end
				end # for: iZ=2:NiZ

		# Averaging	
			if Ngood ≥ 1
				ΔT₂_New = √(ΔT_New_Norm / Float64(Ngood))
			else
				ΔT₂_New = paramHypix.ΔT_Max
			end
			
		return ΔT₂_New, Δθ₂_Max
		end # function ADAPTIVE_TIMESTEP

end # module timeStep
# ...........................................................................................