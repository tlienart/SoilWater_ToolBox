# =============================================================
#		MODULE: bestUniv
# =============================================================
module bestUniv
	import ..option, ..sorptivity, ..wrc, ..kunsat, ..option
	export  BEST_UNIVERSAL_START, TIME_TRANS_STEADY_INDEP, CONVERT_3D_2_1D, CONVERT_1D_2_3D

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BEST_UNIVERSAL_START(iSoil, N_Infilt, ∑Infilt, T, infiltParam, hydroInfilt, infiltOutput)
			# Initializing
				Se_Ini = wrc.θ_2_Se(infiltParam.θ_Ini[iSoil], iSoil, hydroInfilt)

				Kr_θini = (kunsat.Se_2_KUNSAT(Se_Ini, iSoil, hydroInfilt)) / hydroInfilt.Ks[iSoil]

				Sorptivity = sorptivity.SORPTIVITY(infiltParam.θ_Ini[iSoil], iSoil, hydroInfilt)

				A = bestUniv.A(infiltParam.θ_Ini[iSoil], hydroInfilt.θs[iSoil], iSoil, infiltParam)

				B = bestUniv.B(iSoil, Kr_θini, infiltParam)

				T_TransSteady = TIME_TRANS_STEADY(B, hydroInfilt.Ks[iSoil], Sorptivity)

# 				T_TransSteady =   infiltOutput.T_TransSteady_Data[iSoil]

				for iT = 1:N_Infilt[iSoil]
					∑Infilt[iSoil, iT] = BEST_UNIVERSAL(iSoil, A, B, Sorptivity, T[iSoil,iT], T_TransSteady, hydroInfilt, infiltParam)
				end  # for iT=1:N_Infilt[iSoil]

				return ∑Infilt, T_TransSteady
		end # function: BEST_UNIVERSAL_START

		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST_UNIVERSAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BEST_UNIVERSAL(iSoil, A, B, Sorptivity, T, T_TransSteady, hydroInfilt, infiltParam)
			if option.infilt.SingleDoubleRing == "Single" #<>=<>=<>=<>=<>
				if T <= T_TransSteady
					return ∑Infilt = bestUniv.INFILTRATION_3D_TRANSIT(A, B, hydroInfilt.Ks[iSoil], Sorptivity, T)
				else
					return ∑Infilt = bestUniv.INFILTRATION_3D_STEADY(A, B, iSoil, hydroInfilt.Ks[iSoil], Sorptivity, T,  infiltParam, T_TransSteady)
				end # T <= T_TransSteady

			elseif option.infilt.SingleDoubleRing == "Double"  #<>=<>=<>=<>=<>
				if T <= T_TransSteady
					return ∑Infilt = bestUniv.INFILTRATION_1D_TRANSIT(B, hydroInfilt.Ks[iSoil], Sorptivity, T)
				else
					return ∑Infilt = bestUniv.INFILTRATION_1D_STEADY(B, iSoil, hydroInfilt.Ks[iSoil], Sorptivity, T, infiltParam, T_TransSteady)
				end # T <= T_TransSteady
			end # option.∑Infilt.Dimension
		end  # function: BEST_UNIVERSAL



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_3D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_3D_TRANSIT(A, B, Ks, Sorptivity, T)
			return Sorptivity * (T ^ 0.5) + (A * (Sorptivity ^ 2.0) + B * Ks) * T
		end  # function: INFILTRATION_3D_TRANSIT
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_1D_TRANSIT(B, Ks, Sorptivity, T)
			return Sorptivity * (T ^ 0.5) + B * Ks * T
		end # function: INFILTRATION_1D_TRANSIT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_3D_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function INFILTRATION_3D_STEADY(A, B, iSoil, Ks, Sorptivity, T,  infiltParam, T_TransSteady)
		if option.infilt.bestUniv.Continous == true
			return bestUniv.INFILTRATION_3D_TRANSIT(A, B, Ks, Sorptivity, T_TransSteady)  + (A * (Sorptivity ^ 2.0) + Ks) * (T - T_TransSteady)
		elseif option.infilt.bestUniv.Continous  == false
			return (A * (Sorptivity ^ 2.0) + Ks) * T + bestUniv.C(B,  infiltParam, iSoil) * (Sorptivity ^ 2.0) / Ks
		end
	end  # function: INFILTRATION_3D_STEADY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_1D_STEADY(B, iSoil, Ks, Sorptivity, T, infiltParam, T_TransSteady)
			if option.infilt.bestUniv.Continous == true
				return bestUniv.INFILTRATION_1D_TRANSIT(B, Ks, Sorptivity, T_TransSteady)  + Ks * (T - T_TransSteady)
			elseif option.infilt.bestUniv.Continous == false
				return Ks * T + bestUniv.C(B,  infiltParam, iSoil) * (Sorptivity ^ 2.0) / Ks
			end
		end  # function: INFILTRATION_1D_STEADY


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERT_3D_2_1D
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CONVERT_3D_2_1D(hydroInfilt, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)
			for iSoil=1:N_SoilSelect
				
				Sorptivity = sorptivity.SORPTIVITY(infiltParam.θ_Ini[iSoil], iSoil, hydroInfilt)

				A = bestUniv.A(infiltParam.θ_Ini[iSoil], hydroInfilt.θs[iSoil], iSoil, infiltParam)

				for iT=1:N_Infilt[iSoil]
					∑Infilt[iSoil,iT] = ∑Infilt[iSoil,iT] - A * (Sorptivity ^ 2.0) * T[iSoil,iT]
				end # iT
			end # iSoil

			return ∑Infilt
		end  # function: CONVERT_3D_2_1D
			
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERT_3D_2_1D
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CONVERT_1D_2_3D(hydroInfilt, ∑Infilt, infiltParam, N_Infilt, N_SoilSelect, T)
			for iSoil=1:N_SoilSelect

				Sorptivity = sorptivity.SORPTIVITY(infiltParam.θ_Ini[iSoil], iSoil, hydroInfilt)

				A = bestUniv.A(infiltParam.θ_Ini[iSoil], hydroInfilt.θs[iSoil], iSoil, infiltParam)

				for iT=1:N_Infilt[iSoil]
					∑Infilt[iSoil,iT] = ∑Infilt[iSoil,iT] + A * (Sorptivity ^ 2.0) * T[iSoil,iT]
				end # iT
			end # iSoil

			return ∑Infilt
		end  # function: CONVERT_1D_2_3D



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_TRANS_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_TRANS_STEADY(B, Ks, Sorptivity)
			return ( Sorptivity / (Ks * 2.0 * (1.0 - B)) ) ^ 2.0
		end # function: TIME_TRANS_STEADY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : A
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function A(θ_Ini, θs, iSoil, infiltParam)
			return  infiltParam.γ[iSoil] / ( infiltParam.RingRadius[iSoil] * (θs - θ_Ini)) # Units [mm-1]
		end  # function: A
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : B
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function B(iSoil, Kr_θini,  infiltParam)
			return (2.0 -  infiltParam.β[iSoil]) / 3.0 + Kr_θini * (1.0 + infiltParam.β[iSoil]) / 3.0
		end # function: B


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : C
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function C(B, infiltParam, iSoil)
			return log(1.0 /  infiltParam.β[iSoil]) * (1.0 +  infiltParam.β[iSoil]) / (6.0 * (1.0 -  infiltParam.β[iSoil]) * (1.0 - B) )
		end # function: C

end  # module: bestUniv
# ............................................................