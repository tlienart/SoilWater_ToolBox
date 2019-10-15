# =============================================================
#		MODULE: best
# =============================================================
module best
	using ..option, ..sorptivity, ..wrc, ..kunsat
	export BEST

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BEST(iSoil, Dimension, Sorptivity, Tinfilt, θ_Ini, hydro, infilt)

			# Required data
				Se_Ini = wrc.θ_2_Se(θ_Ini, iSoil, hydro)

				Kr_θini= (kunsat.Se_2_KUNSAT(Se_Ini, iSoil, hydro)) / hydro.Ks[iSoil]

			# Best parameters
				A = best.A(θ_Ini, hydro.θs[iSoil], iSoil, infilt)

				B = best.B(iSoil, Kr_θini, infilt)

				Time_TransStead = TIME_TRANS_STEADY(B, hydro.Ks[iSoil], Sorptivity)

			# Required options 1D or 3D
				if Dimension	== "3D"
					if Tinfilt <= Time_TransStead 
						return Infilt = best.INFILTRATION_3D_TRANSIT(A, B, hydro.Ks[iSoil], Sorptivity, Tinfilt)
					else
						return Infilt = best.INFILTRATION_3D_STEADY(A, B, iSoil, hydro.Ks[iSoil], Sorptivity, Tinfilt, infilt)
					end # Tinfilt <= Time_TransStead 

				elseif Dimension == "1D"
					if Tinfilt <= Time_TransStead 
						return Infilt = best.INFILTRATION_1D_TRANSIT(B, hydro.Ks[iSoil], Sorptivity, Tinfilt)
					else
						return Infilt = best.INFILTRATION_1D_STEADY(B, iSoil, hydro.Ks[iSoil], Sorptivity, Tinfilt, infilt)
					end # Tinfilt <= Time_TransStead 
				end # option.Infilt.Dimension
		end # function: BEST


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_3D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_3D_TRANSIT(A, B, Ks, Sorptivity, Tinfilt)
			return Sorptivity * (Tinfilt ^ 0.5) + (A * (Sorptivity ^ 2.0) + B * Ks) * Tinfilt
		end  # function: INFILTRATION_3D_TRANSIT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_3D_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_3D_STEADY(A, B, iSoil, Ks, Sorptivity, Tinfilt, infilt)
			return (A * (Sorptivity ^ 2.0) + Ks) * Tinfilt + C(B, infilt, iSoil) * (Sorptivity ^ 2.0) / Ks
		end  # function:


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_1D_TRANSIT(B, Ks, Sorptivity, Tinfilt)
			return Sorptivity * (Tinfilt ^ 0.5) + B * Ks * Tinfilt
		end # function: INFILTRATION_1D_TRANSIT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_1D_STEADY(B, iSoil, Ks, Sorptivity, Tinfilt, infilt)
			return Ks * Tinfilt + best.C(B, infilt, iSoil) * (Sorptivity^2.) / Ks
		end  # function: INFILTRATION_1D_STEADY


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_TRANS_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_TRANS_STEADY(B, Ks, Sorptivity)
			return ( Sorptivity / (Ks * 2.0 * (1.0 - B)) ) ^ 2.0
		end # function: TIME_TRANS_STEADY

		
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_TRANS_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_TRANS_STEADY_INDEP(iSoil, θ_Ini, hydro, infilt)
			Kr_θini = (kunsat.θ_2_KUNSAT(θ_Ini, iSoil, hydro)) / hydro.Ks[iSoil]
			
			B = best.B(iSoil, Kr_θini, infilt)

			Sorptivity = sorptivity.SORPTIVITY(θ_Ini, iSoil, hydro)

			return ( Sorptivity / (hydro.Ks[iSoil] * 2.0 * (1.0 - B)) ) ^ 2.0
		end # function: TIME_TRANS_STEADY


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : A
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function A(θ_Ini, θs, iSoil, infilt)
			return infilt.γ[iSoil] / ( infilt.RingRadius[iSoil] * (θs - θ_Ini) ) # Units (mm-1)
		end  # function: A
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : B
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function B(iSoil, Kr_θini, infilt)
			return (2.0 - infilt.β[iSoil]) / 3.0 + Kr_θini * (1.0 + infilt.β[iSoil]) / 3.0
		end # function: B


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : C
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function C(B, infilt, iSoil)
			return log(1.0 / infilt.β[iSoil]) * (1.0 + infilt.β[iSoil]) / (6.0 * (1.0 - infilt.β[iSoil]) * (1.0 - B) )
		end # function: C

end  # module: best
# ............................................................