module interpolate
	export ∑_2_Δ, POINTS_2_SlopeIntercept, INTERPOLATE_2D_LOOP

	"""
	∑RAIN_2_ΔPr(∑Pr, ∑Pr_Climate, ∑T_Climate, ∑T, iT, iT_Pr, N_Climate; FlagForwardTime=true)
	This is used for other variables, we give a example for Pr
	"""
	function ∑_2_Δ(∑X_Past, ∑X_Climate, ∑T, ∑T_Climate, iT_X, N_Climate, Flag_ReRun, iT)

		# Moving backwards if we need to rerun
			if Flag_ReRun && iT_X ≥ 3
				iT_X -= 1
			end

		# Determening if we should increase iT_X
			FlagBreak = false
			while !(FlagBreak)
				if (∑T_Climate[iT_X-1] ≤ ∑T[iT] ≤ ∑T_Climate[iT_X]) || (iT_X == N_Climate) 
					FlagBreak = true
					break
				else 
					iT_X += 1
					FlagBreak = false
				end # if
			end # while

		# Building a regression line which passes from POINT1(∑T_Climate[iT_X], ∑Pr_Climate[iT_Pr]) and POINT2: (∑T_Climate[iT_Pr+1], ∑Pr_Climate[iT_Pr+1])
			Slope, Intercept = POINTS_2_SlopeIntercept(∑T_Climate[iT_X-1], ∑X_Climate[iT_X-1], ∑T_Climate[iT_X], ∑X_Climate[iT_X])

			∑X = Slope * ∑T[iT] + Intercept

		# Xecipitation [mm /  ΔTconst]
			ΔX = ∑X - ∑X_Past
		
	return ∑X, ΔX, iT_X
	end # function ∑RAIN_2_ΔX


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : POINTS_2_SlopeIntercept
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	"""POINTS_2_SlopeIntercept
	From Point1 [X1, Y1] and point2 [X2, Y2] compute Y = Slope.X₀ + Intercept
	"""
		function POINTS_2_SlopeIntercept(X1, Y1, X2, Y2)
			Slope = (Y2 - Y1) / (X2 - X1)
			Intercept = (Y1 * X2 - X1 * Y2) / (X2 - X1)
			return Slope, Intercept
		end # POINTS_2_SlopeIntercept


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INTERPOLATE_2D_LOOP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INTERPOLATE_2D_LOOP(∑T, ∑T_Reduced, Nit, NiZ, X₀_Reduced, X₀)

			Nit_Reduced = length(∑T_Reduced)
			iT_X = 2

			for iT_Reduced = 1:Nit_Reduced
				FlagBreak = false

				while !(FlagBreak)
					if ∑T_Reduced[iT_Reduced] < ∑T[iT_X-1]
						error("HYPIX_MODEL INTERPOLATE_2D_LOOP:  ∑T_Reduced[iT_Reduced] < ∑T[iT_X-1] iT_Reduced=$iT_Reduced iT_X=$iT_X")
					end

					if (∑T[iT_X-1] ≤ ∑T_Reduced[iT_Reduced] ≤ ∑T[iT_X]) || (iT_X == Nit) 
						FlagBreak = true
						break

					else 
						iT_X += 1
						FlagBreak = false
					end # if
				end # while

				# Building a regression line which passes from POINT1(∑T_Climate[iT_X], ∑Pr_Climate[iT_Pr]) and POINT2: (∑T_Climate[iT_Pr+1], ∑Pr_Climate[iT_Pr+1])
				# if !Flag_TooEarly
					for iZ = 1:NiZ
						Slope, Intercept = interpolate.POINTS_2_SlopeIntercept(∑T[iT_X-1], X₀[iT_X-1,iZ], ∑T[iT_X], X₀[iT_X,iZ])

						X₀_Reduced[iT_Reduced,iZ] = Slope * ∑T_Reduced[iT_Reduced] + Intercept
					end # for iZ = 1:NiZ	
			end # for: iT_Reduced=1:obsTheta.Nit
				
	return X₀_Reduced
	end  # function: θINTERPOLATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INTERPOLATE_2D_LOOP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, X₀_Reduced, X₀)
			iT_X = 2
			for iT_Reduced=1:Nit_Reduced
		
				FlagBreak = false
				while !(FlagBreak)
					if (∑T[iT_X-1] - eps(10.0) ≤ ∑T_Reduced[iT_Reduced] ≤ ∑T[iT_X] + eps(10.0)) || (iT_X == Nit) 
						FlagBreak = true
						break
					else 
						iT_X += 1
						FlagBreak = false
					end # if
				end # while

				# Building a regression line which passes from POINT1(∑T_Climate[iT_X], ∑Pr_Climate[iT_Pr]) and POINT2: (∑T_Climate[iT_Pr+1], ∑Pr_Climate[iT_Pr+1])
				Slope, Intercept = interpolate.POINTS_2_SlopeIntercept(∑T[iT_X-1], X₀[iT_X-1], ∑T[iT_X], X₀[iT_X])

				X₀_Reduced[iT_Reduced] = Slope * ∑T_Reduced[iT_Reduced] + Intercept
			
			end # for: iT_Reduced=1:obsTheta.Nit
		
		return X₀_Reduced
	end  # function: θINTERPOLATION


#	 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INTERPOLATE_1D_MAX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function INTERPOLATE_1D_MAX(∑T, ∑T_Reduced, Nit_Reduced, Nit::Int64, X₀_Reduced, X₀)
		iT_X = 2::Int64
		for iT_Reduced=1:Nit_Reduced
			FlagBreak = false
			Xmax = 0.0::Float64
			while !(FlagBreak)
				Xmax = max(Xmax,  X₀[iT_X])
				if (∑T[iT_X-1] - eps(10.0) ≤ ∑T_Reduced[iT_Reduced] ≤ ∑T[iT_X] + eps(10.0)) || (iT_X == Nit) 
					FlagBreak = true
					break
				else 
					iT_X += 1
					FlagBreak = false
				end # if
			end # while

			X₀_Reduced[iT_Reduced] = Xmax
		end # for: iT_Reduced=1:obsTheta.Nit
	
	return X₀_Reduced
	end  # function: INTERPOLATE_1D_MAX

end # module interpolate