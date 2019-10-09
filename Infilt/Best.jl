# =============================================================
#		MODULE: best
# =============================================================
module best
	include ..option


	function BEST(Time, θ_Ini, iSoil, Sorptivity, Kr_θini, Ks, hydro, infilt, )
		B = B(Kr_θini, infilt, iSoil)

		A = A(θ_Ini, iSoil, hydro, infilt)

		Time_TransStead =  TIME_TRANS_STEADY(Sorptivity, B, Ks)

		if Time <= Time_TransStead 
			return Infilt = INFILTRATION_3D_TRANSIT(Time, Sorptivity, A, B, Ks)
		else
			if Flag_Best == "Best_G"
				return Infilt = INFILTRATION_3D_STEADY_BESTG(Time, Sorptivity, A, B, Ks)
			else
				return Infilt = INFILTRATION_3D_STEADY_BESTGI(Time, Time_TransStead, Sorptivity, A, B, Ks)
			end
		end
	end

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_3D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_3D_TRANSIT(Time, Sorptivity, A, B, Ks)
			return Sorptivity * (Time ^ 0.5) + (A * (Sorptivity ^ 2.0) + B * Ks) * Time
		end  # function: INFILTRATION_3D_TRANSIT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_3D_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_3D_STEADY(Time, Sorptivity, A, B, Ks)
			return (A * (Sorptivity ^ 2.0) + Ks) * Time + C(B) * (Sorptivity ^ 2.0) / Ks
		end  # function:


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : A
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function A(θ_Ini, iSoil, hydro, infilt)
			return infilt.γ[iSoil] / ( infilt.RingRadius[iSoil] * (hydro.θs[iSoil] - θ_Ini) ) # Units (mm-1)
		end  # function: A
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : B
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function B(Kr_θini, infilt, iSoil)
			return (2.0 - infilt.β[iSoil]) / 3.0 + Kr_θini * (1.0 + infilt.β[iSoil]) / 3.0
		end # function: B


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : C
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function C(B, infilt, iSoil)
			return log(1.0 / infilt.β[iSoil]) * (1.0 + infilt.β[iSoil]) / (6. * (1.0 - infilt.β[iSoil]) * (1.0 - B) )
		end # function: C


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_TRANS_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_TRANS_STEADY(Sorptivity, B, Ks)
			return ( Sorptivity / (Ks * 2.0 * (1.0 - B)) ) ^ 2.0
		end # function: TIME_TRANS_STEADY


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_1D_TRANSIT(Time, Sorptivity, B, Ks)
			return Sorptivity * (Time ^ 0.5) + B * Ks * Time
		end # function: INFILTRATION_1D_TRANSIT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_1D_STEADY(Time, Sorptivity, A, B, Ks)
			return Ks * Time + C(B) * (Sorptivity^2.) / Ks
		end  # function: INFILTRATION_1D_STEADY




	function Q1D_TRANSIT(Time, Sorptivity, Se_Ini, θs, θr, σ, Ks)
		Kr_θini = kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, 1., θs, σ)
		Q1D = 0.5 * Sorptivity * (Time^-0.5) + B(Kr_θini)*Ks
		return Q1D
	end








end  # module: best
# ............................................................