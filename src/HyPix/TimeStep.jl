module timeStep
	import ..wrc
   export TIMESTEP, ADAPTIVE_TIMESTEP, ΔΨMAX!

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION :  TIMESTEP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIMESTEP(∑T, discret, Flag_ReRun::Bool, hydro, iT::Int64, N_∑T_Climate::Float64, NiZ::Int64, option, param, Q, ΔLnΨmax, ΔSink, ΔT, θ, Ψ)

			Δθ_Max = param.hyPix.Δθ_Max

			# The iT is of the previous simulation
			if !Flag_ReRun # <>=<>=<>=<>=<>	
				ΔT₂, Δθ_Max = ADAPTIVE_TIMESTEP(discret, hydro, iT, NiZ, option, param, Q, ΔLnΨmax, ΔSink, θ, Ψ)
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
		function ΔΨMAX!(hydro, NiZ::Int64, option, param, ΔLnΨmax::Vector{Float64})

			for iZ=1:NiZ
				θ½ = (hydro.θsMacMat[iZ] + hydro.θr[iZ]) * 0.5
				
				θ△ = min(θ½ + param.hyPix.Δθ_Max * 0.5, hydro.θs[iZ])

				θ▽ = max(θ½ - param.hyPix.Δθ_Max * 0.5, hydro.θr[iZ])

				ΔLnΨmax[iZ] = (log(wrc.θ_2_ΨDual(option.hyPix, θ▽, iZ, hydro)) - log(wrc.θ_2_ΨDual(option.hyPix, θ△, iZ, hydro))) * 0.5			
			end # for iZ=1:NiZ
	
		return ΔLnΨmax
		end  # function: ΔΨMAX!
	#--------------------------------------------------------------------


	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΔθMAX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ΔθMAX(hydro, iT::Int64, iZ::Int64, option, ΔLnΨmax::Vector{Float64}, Ψ::Matrix{Float64})
			
			if Ψ[iT,iZ] > ΔLnΨmax[iZ]
				Ψ▽ = exp(log(Ψ[iT,iZ]) - ΔLnΨmax[iZ])
			else
				Ψ▽ = 0.0
			end

			Ψ△  = exp(log(Ψ[iT,iZ]) + ΔLnΨmax[iZ])
			
			θ△ = wrc.Ψ_2_θDual(option.hyPix, Ψ▽, iZ, hydro)
		
			θ▽ = wrc.Ψ_2_θDual(option.hyPix, Ψ△, iZ, hydro)
	
		return θ△ - θ▽
		end  # function:  ΔθMAX
	# ------------------------------------------------------------------



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ADAPTIVE_TIMESTEP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ADAPTIVE_TIMESTEP(discret, hydro, iT::Int64, NiZ::Int64, option, param, Q, ΔLnΨmax, ΔSink, θ, Ψ)

			# Searching for the minimum value of ΔT of the simulation
				if option.hyPix.NormMin⍰ == "Norm"
					ΔT_New_Norm = 0.0
				else
					ΔT_New_Norm = Inf
				end
			
			# Initializing
				Δθ₂_Max = param.hyPix.Δθ_Max

			# Computing smallest Δθ_Max
				Ngood = 0
				for iZ = 1:NiZ-1

					if abs(Ψ[iT,iZ] - Ψ[iT,iZ+1]) ≥ 1.0

						if option.hyPix.AdaptiveTimeStep⍰ == "ΔΨ" # <>=<>=<>=<>=<>
							Δθ₂_Max = ΔθMAX(hydro, iT, iZ, option, ΔLnΨmax, Ψ)	
						end # option.hyPix.AdaptiveTimeStep⍰ ==:ΔΨ
					
						ΔT₂_New = (discret.ΔZ[iZ] * Δθ₂_Max + ΔSink[iT,iZ]) / (abs(Q[iT,iZ] - Q[iT,iZ+1]))

						ΔT₂_New = min( max(param.hyPix.ΔT_Min, ΔT₂_New), param.hyPix.ΔT_Max)
		
						if option.hyPix.NormMin⍰ == "Norm"
							ΔT_New_Norm += ΔT₂_New ^ 2
						else
							ΔT_New_Norm = min(ΔT_New_Norm, ΔT₂_New)
						end
						Ngood += 1
					end

				end # for: iZ=2:NiZ

		# Averaging	
			if option.hyPix.NormMin⍰ == "Norm"
				if Ngood ≥ 1
					ΔT₂_New = √(ΔT_New_Norm / Float64(Ngood))
				else
					ΔT₂_New = param.hyPix.ΔT_Max
				end
			else
				if Ngood ≥ 1
					ΔT₂_New = ΔT_New_Norm
				else
					ΔT₂_New = param.hyPix.ΔT_Max
				end
			end	
		return ΔT₂_New, Δθ₂_Max
		end # function ADAPTIVE_TIMESTEP

end # module timeStep
# ...........................................................................................