# =============================================================
#		module: kunsatModel
# =============================================================
module θψ2KsModel
	import ..cst, ..distribution
	import QuadGK
	import SpecialFunctions: erfc, erfcinv
	export KSMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KS_MODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSMODEL(hydro, KₛModel, KₛModel⍰, ksmodelτ, N_iZ, optim, optimKsmodel; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])
			for iZ=1:N_iZ
				if Flag_RockFragment
					RockFragment₁ = RockFragment[iZ]
				else
					RockFragment₁ = 0.0
				end #@isdefined RockFragment

				if Flag_IsTopsoil
					IsTopsoil₁ = Int64(IsTopsoil[iZ])
				else
					IsTopsoil₁ = 1	# Default value				
				end  # if: @isdefined IsTopsoil

				KₛModel[iZ] = θΨ_2_KSMODEL(hydro, IsTopsoil₁, iZ,  KₛModel⍰, ksmodelτ; RockFragment=RockFragment₁)
			end # if: hydro.Ks[iZ] > eps(10.0)	
		return KₛModel
		end  # function: KS_MODEL
	#..................................................................

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : KUNSAT_MODEL
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""Pollacco, J.A.P., Webb, T., McNeill, S., Hu, W., Carrick, S., Hewitt, A., Lilburne, L., 2017. Saturated hydraulic conductivity model computed from bimodal water retention curves for a range of New Zealand soils. Hydrol. Earth Syst. Sci. 21, 2725–2737. https://doi.org/10.5194/hess-21-2725-2017"""
		function θΨ_2_KSMODEL(hydroParam, IsTopsoil, iZ::Int64, KₛModel⍰,ksmodelτ; RockFragment=0.0, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ], τ₁=ksmodelτ.τ₁[IsTopsoil], τ₂=ksmodelτ.τ₂[IsTopsoil], τ₃=ksmodelτ.τ₃[IsTopsoil], τ₄=ksmodelτ.τ₄[IsTopsoil], τ₅=ksmodelτ.τ₅[IsTopsoil], τ₁Mac=ksmodelτ.τ₁Mac[IsTopsoil], τ₂Mac=ksmodelτ.τ₂Mac[IsTopsoil], τ₃Mac=ksmodelτ.τ₃Mac[IsTopsoil], RockFragment_Treshold=0.4, Rtol=1.0E-3, Se_Max=0.9999 )

			# Determine when Ks increases for increasing RockFragment	
				if RockFragment > RockFragment_Treshold

					RockFragment2 = max(2.0 * RockFragment_Treshold - RockFragment, 0.0)

					θs = (θs / (1.0 - RockFragment)) * (1.0 - RockFragment2)
					
					θsMacMat = (θsMacMat / (1.0 - RockFragment)) * (1.0 - RockFragment2)

					θr = (θr / (1.0 - RockFragment)) * (1.0 - RockFragment2)
				end				

			# Model traditional 1 =================================================================
			if  KₛModel⍰=="Model1" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				# Transformation matrix
					T1 =  (10.0 ^ - (1.0 / (1.0 - τ₁ ))) / 0.1
					T2 = 2.0 * (1.0 - τ₂)
					T3 = 1.0 / (1.0 - τ₃)

				# Transformation macro
					T1Mac = (10.0 ^ - (1.0 / (1.0 - τ₁ *τ₁Mac))) / 0.1
					T2Mac = 2.0 * (1.0 - τ₂Mac * τ₂)
					T3Mac = 1.0 / (1.0 - τ₃Mac * τ₃)

					# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					#		FUNCTION : KS_MODEL_1
					# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
						function KS_MODEL_1(Se, T1, T2, T3, T1Mac, T2Mac, T3Mac)

							Kunsat_Matrix =  T1 * ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

							if θs - θsMacMat > 0.001
								Kunsat_Macro = T1Mac * ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
							else
								Kunsat_Macro = 0.0
							end
						return Kunsat = Kunsat_Matrix + Kunsat_Macro
						end  # function: KS_MODEL_1
					# ------------------------------------------------------------------
				
					KₛModel = cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_1(Se, T1, T2, T3, T1Mac, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]

				return KₛModel

			elseif  KₛModel⍰=="Model2"  # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				# Transformation matrix
					T1 =  (10.0 ^ - (1.0 / (1.0 - τ₁ ))) / 0.1
					T2 = 2.0 * (1.0 - τ₂)
					T3 = 1.0 / (1.0 - τ₃)

				# Transformation macro
					T2Mac = 2.0 * (1.0 - τ₂Mac * τ₂)
					T3Mac = 1.0 / (1.0 - τ₃Mac * τ₃)

				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				#		FUNCTION : KS_MODEL_2
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					function KS_MODEL_2(Se, T2, T3, T2Mac, T3Mac)
						Kunsat_Matrix =   ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

						if θs - θsMacMat > 0.001
							Kunsat_Macro = ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
						else
							Kunsat_Macro = 0.0
						end
					return Kunsat = Kunsat_Matrix + Kunsat_Macro
					end  # function: KS_MODEL_1
				# ------------------------------------------------------------------
			
				KₛModel = T1 * cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_2(Se, T2, T3, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]
				return KₛModel
				
			elseif  KₛModel⍰=="Model3"  # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>

				# Transformation matrix
					Tσ = τMODEL_σ(hydroParam, iZ, τ₁, τ₁ * τ₄, σ)

					T1 =  (10.0 ^ - (1.0 / (1.0 - Tσ ))) / 0.1
					T2 = 2.0 * (1.0 - τ₂)
					T3 = 1.0 / (1.0 - τ₃)

				# Transformation macro
					T2Mac = 2.0 * (1.0 - τ₂Mac * τ₂)
					T3Mac = 1.0 / (1.0 - τ₃Mac * τ₃)

				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				#		FUNCTION : KS_MODEL_1
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function KS_MODEL_3(Se, T2, T3, T2Mac, T3Mac)

					Kunsat_Matrix =  ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

					if θs - θsMacMat > 0.001
						Kunsat_Macro = ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
					else
						Kunsat_Macro = 0.0
					end
				return Kunsat = Kunsat_Matrix + Kunsat_Macro
				end  # function: KS_MODEL_3
			# ------------------------------------------------------------------
			
				KₛModel =  T1 * cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_3(Se, T2, T3, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]
				return KₛModel
			
			
			elseif  KₛModel⍰=="Model4"  # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				# Transformation matrix

				Tσ = τMODEL_σSilt(hydroParam, iZ, τ₁, τ₁ * τ₄, σ; σSilt_η=0.538, Pσ=3.0, Distribution⍰="Normal", Normalise=true, Invert=false)

				T1 =  (10.0 ^ - (1.0 / (1.0 - Tσ ))) / 0.1
				T2 = 2.0 * (1.0 - τ₂)
				T3 = 1.0 / (1.0 - τ₃)

			# Transformation macro
				T2Mac = 2.0 * (1.0 - τ₂Mac * τ₂)
				T3Mac = 1.0 / (1.0 - τ₃Mac * τ₃)

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : KS_MODEL_1
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function KS_MODEL_4(Se, T2, T3, T2Mac, T3Mac)

				Kunsat_Matrix =  ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

				if θs - θsMacMat > 0.001
					Kunsat_Macro = ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
				else
					Kunsat_Macro = 0.0
				end
			return Kunsat = Kunsat_Matrix + Kunsat_Macro
			end  # function: KS_MODEL_4
		# ------------------------------------------------------------------
		
			KₛModel =  T1 * cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_4(Se, T2, T3, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]
			return KₛModel
				
		elseif  KₛModel⍰=="Model5"  # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
			# Transformation matrix

			Tσ = τMODEL_σSilt_σ(hydroParam, iZ, τ₁,  τ₁ * τ₄, σ; Amplitude=τ₅, σSilt_η=0.538, Pσ=3, Distribution⍰="Normal", Normalise=true, Invert=false)

			T1 =  (10.0 ^ - (1.0 / (1.0 - Tσ ))) / 0.1
			T2 = 2.0 * (1.0 - τ₂)
			T3 = 1.0 / (1.0 - τ₃)

		# Transformation macro
			T2Mac = 2.0 * (1.0 - τ₂Mac * τ₂)
			T3Mac = 1.0 / (1.0 - τ₃Mac * τ₃)

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KS_MODEL_1
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KS_MODEL_5(Se, T2, T3, T2Mac, T3Mac)

			Kunsat_Matrix =  ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

			if θs - θsMacMat > 0.001
				Kunsat_Macro = ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
			else
				Kunsat_Macro = 0.0
			end
		return Kunsat = Kunsat_Matrix + Kunsat_Macro
		end  # function: KS_MODEL_4
	# ------------------------------------------------------------------
	
		KₛModel =  T1 * cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_5(Se, T2, T3, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]
		return KₛModel
			
		elseif  KₛModel⍰=="Model6"  # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>

			Tσ = τMODEL_σSilt(hydroParam, iZ, τ₂, τ₂ * τ₄, σ; σSilt_η=τ₅, Pσ=3.0, Distribution⍰="Normal", Normalise=true, Invert=false)

			T1 =  (10.0 ^ - (1.0 / (1.0 - τ₁ ))) / 0.1
			T2 = 2.0 * (1.0 - Tσ)
			T3 = 1.0 / (1.0 - τ₃)

		# Transformation macro
			T2Mac = 2.0 * (1.0 - τ₂Mac * Tσ * τ₂)
			T3Mac = 1.0 / (1.0 - τ₃Mac * τ₃)

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KS_MODEL_1
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KS_MODEL_6(Se, T2, T3, T2Mac, T3Mac)

			Kunsat_Matrix =  ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

			if θs - θsMacMat > 0.001
				Kunsat_Macro = ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
			else
				Kunsat_Macro = 0.0
			end
		return Kunsat = Kunsat_Matrix + Kunsat_Macro
		end  # function: KS_MODEL_6
	# ------------------------------------------------------------------
	
		KₛModel =  T1 * cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_6(Se, T2, T3, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]
		return KₛModel
		

		elseif  KₛModel⍰=="Model7"  #sre <>=<>=<>=<>=<>=<>=<>=<>=<>=<>

			Tσ = τMODEL_σ(hydroParam, iZ, τ₂, τ₂ * τ₄, σ)

			T1 =  (10.0 ^ - (1.0 / (1.0 - τ₁ ))) / 0.1
			T2 = 2.0 * (1.0 - Tσ)
			T3 = 1.0 / (1.0 - τ₃)

		# Transformation macro
			T2Mac = 2.0 * (1.0 - τ₂Mac * τ₂)
			T3Mac = 1.0 / (1.0 - τ₃Mac * τ₃)

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KS_MODEL_1
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KS_MODEL_7(Se, T2, T3, T2Mac, T3Mac)

			Kunsat_Matrix =  ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

			if θs - θsMacMat > 0.001
				Kunsat_Macro = ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
			else
				Kunsat_Macro = 0.0
			end
		return Kunsat = Kunsat_Matrix + Kunsat_Macro
		end  # function: KS_MODEL_6
	# ------------------------------------------------------------------
	
		KₛModel =  T1 * cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_7(Se, T2, T3, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]
		return KₛModel

		elseif  KₛModel⍰=="Model8"  #sre <>=<>=<>=<>=<>=<>=<>=<>=<>=<>

			Tσ = τMODEL_σ(hydroParam, iZ, τ₃, τ₃ * τ₄, σ)

			T1 =  (10.0 ^ - (1.0 / (1.0 - τ₁ ))) / 0.1
			T2 = 2.0 * (1.0 - τ₂)
			T3 = 1.0 / (1.0 - Tσ)

		# Transformation macro
			T2Mac = 2.0 * (1.0 - τ₂Mac * τ₂)
			T3Mac = 1.0 / (1.0 - τ₃Mac * Tσ)

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KS_MODEL_1
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KS_MODEL_8(Se, T2, T3, T2Mac, T3Mac)

			Kunsat_Matrix =  ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

			if θs - θsMacMat > 0.001
				Kunsat_Macro = ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
			else
				Kunsat_Macro = 0.0
			end
		return Kunsat = Kunsat_Matrix + Kunsat_Macro
		end  # function: KS_MODEL_6
	# ------------------------------------------------------------------
	
		KₛModel =  T1 * cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_8(Se, T2, T3, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]
		return KₛModel

		elseif  KₛModel⍰=="Model9"  #sre <>=<>=<>=<>=<>=<>=<>=<>=<>=<>

			Tσ = τMODEL_σSilt(hydroParam, iZ, τ₃, τ₃ * τ₄, σ; σSilt_η=0.538, Pσ=3.0, Distribution⍰="Normal", Normalise=true, Invert=false)
	
			T1 =  (10.0 ^ - (1.0 / (1.0 - τ₁ ))) / 0.1
			T2 = 2.0 * (1.0 - τ₂)
			T3 = 1.0 / (1.0 - Tσ)

		# Transformation macro
			T2Mac = 2.0 * (1.0 - τ₂Mac * τ₂)
			T3Mac = 1.0 / (1.0 - τ₃Mac * Tσ)

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KS_MODEL_1
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KS_MODEL_9(Se, T2, T3, T2Mac, T3Mac)

			Kunsat_Matrix =  ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

			if θs - θsMacMat > 0.001
				Kunsat_Macro = ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
			else
				Kunsat_Macro = 0.0
			end
		return Kunsat = Kunsat_Matrix + Kunsat_Macro
		end  # function: KS_MODEL_6
	# ------------------------------------------------------------------
	
		KₛModel =  T1 * cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_9(Se, T2, T3, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]
		return KₛModel


	elseif  KₛModel⍰=="Model10"  #sre <>=<>=<>=<>=<>=<>=<>=<>=<>=<>

		Tσ1 = τMODEL_σSilt(hydroParam, iZ, τ₁, τ₁ * τ₄, σ; σSilt_η=0.538, Pσ=4.0, Distribution⍰="Normal", Normalise=true, Invert=false)

		Tσ2 = τMODEL_σSilt(hydroParam, iZ, τ₂, τ₂ * τ₅, σ; σSilt_η=0.538, Pσ=4.0, Distribution⍰="Normal", Normalise=true, Invert=false)

		T1 =  (10.0 ^ - (1.0 / (1.0 - Tσ1 ))) / 0.1
		T2 = 2.0 * (1.0 - Tσ2)
		T3 = 1.0 / (1.0 - τ₃)

	# Transformation macro
		T2Mac = 2.0 * (1.0 - τ₂Mac * Tσ2)
		T3Mac = 1.0 / (1.0 - τ₃Mac * τ₃)

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KS_MODEL_1
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function KS_MODEL_10(Se, T2, T3, T2Mac, T3Mac)

		Kunsat_Matrix =  ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

		if θs - θsMacMat > 0.001
			Kunsat_Macro = ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
		else
			Kunsat_Macro = 0.0
		end
	return Kunsat = Kunsat_Matrix + Kunsat_Macro
	end  # function: KS_MODEL_10
# ------------------------------------------------------------------

	KₛModel =  T1 * cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_10(Se, T2, T3, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]
	return KₛModel


elseif  KₛModel⍰=="Model11"  #sre <>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	Tσ1 = τMODEL_σ(hydroParam, iZ, τ₁, τ₁ * τ₄, σ)
	Tσ2 = τMODEL_σ(hydroParam, iZ, τ₂, τ₂ * τ₅, σ)

	T1 =  (10.0 ^ - (1.0 / (1.0 - Tσ1 ))) / 0.1
	T2 = 2.0 * (1.0 - Tσ2)
	T3 = 1.0 / (1.0 - τ₃)

# Transformation macro
	T2Mac = 2.0 * (1.0 - τ₂Mac * Tσ2)
	T3Mac = 1.0 / (1.0 - τ₃Mac * τ₃)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : KS_MODEL_1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function KS_MODEL_11(Se, T2, T3, T2Mac, T3Mac)

	Kunsat_Matrix =  ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

	if θs - θsMacMat > 0.001
		Kunsat_Macro = ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
	else
		Kunsat_Macro = 0.0
	end
return Kunsat = Kunsat_Matrix + Kunsat_Macro
end  # function: KS_MODEL_10
# ------------------------------------------------------------------

KₛModel =  T1 * cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_11(Se, T2, T3, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]
return KₛModel


elseif  KₛModel⍰=="Model12"  #sre <>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	Tσ1 = τMODEL_σ(hydroParam, iZ, τ₁, τ₁ * τ₄, σ)
	Tσ2 = τMODEL_σSilt(hydroParam, iZ, τ₂, τ₂ * τ₅, σ; σSilt_η=0.538, Pσ=2.0, Distribution⍰="Normal", Normalise=true, Invert=false)

	T1 =  (10.0 ^ - (1.0 / (1.0 - Tσ1 ))) / 0.1
	T2 = 2.0 * (1.0 - Tσ2)
	T3 = 1.0 / (1.0 - τ₃)

# Transformation macro
	T2Mac = 2.0 * (1.0 - τ₂Mac * Tσ2)
	T3Mac = 1.0 / (1.0 - τ₃Mac * τ₃)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : KS_MODEL_1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function KS_MODEL_12(Se, T2, T3, T2Mac, T3Mac)

	Kunsat_Matrix =  ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

	if θs - θsMacMat > 0.001
		Kunsat_Macro = ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
	else
		Kunsat_Macro = 0.0
	end
return Kunsat = Kunsat_Matrix + Kunsat_Macro
end  # function: KS_MODEL_10
# ------------------------------------------------------------------

KₛModel =  T1 * cst.KunsatModel * QuadGK.quadgk(Se -> KS_MODEL_12(Se, T2, T3, T2Mac, T3Mac), 0.0, Se_Max; rtol=Rtol)[1]
return KₛModel



		end  # if: Model=="Model3"
	end  # function: θΨ_2_KSMODEL 


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  τMODEL_σSilt
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_σSilt_σ(hydroParam, iZ, Pₘₐₓ, Pₘᵢₙ, σ; Amplitude=0.5, σSilt_η=0.538, Pσ=3, Distribution⍰="Normal", Normalise=true, Invert=false)

			ση = τMODEL_σ(hydroParam, iZ, Pₘₐₓ, Pₘᵢₙ, σ)

			τσ_Dist = τMODEL_σSilt(hydroParam, iZ, Pₘₐₓ, Pₘᵢₙ, ση; σSilt_η=σSilt_η, Pσ=Pσ, Distribution⍰="Normal", Normalise=Normalise, Invert=Invert)

			τσ = τMODEL_σ(hydroParam, iZ, Pₘₐₓ, Pₘᵢₙ, σ)

			return τ = min(τσ + Amplitude * (τσ_Dist / (σSilt_η + 1.0)) , 1.0) * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
		end
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  τMODEL_σSilt
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_σSilt(hydroParam, iZ, Pₘₐₓ, Pₘᵢₙ, σ; σSilt_η=0.538, Pσ=3.0, Distribution⍰="Normal", Normalise=true, Invert=false)

			ση = τMODEL_σ(hydroParam, iZ, Pₘₐₓ, Pₘᵢₙ, σ)

			if  Distribution⍰== "Normal"
				σ_Dist = σSilt_η / Pσ

			elseif  Distribution⍰== "LogNormal"
				σ_Dist = log(σSilt_η) / Pσ

			else
				error("*** τMODEL_σSilt: $Distribution⍰ not implemented try <Normal> or  <LogNormal>  ***")
			end

			τσ_Dist = distribution.DISTRIBUTION(ση, σSilt_η, σ_Dist; Distribution⍰=Distribution⍰, Normalise=Normalise, Invert=Invert)[1]
			return τ = τσ_Dist  * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
		end
		#..................................................................


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : τMODEL_MIXING
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function τMODEL_σ(hydroParam, iZ, Pₘₐₓ, Pₘᵢₙ, σ)
				ση = σ_2_ση(hydroParam, iZ, σ)
				return τ = ση  * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
			end
		#..................................................................


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : σ_2_ση
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function σ_2_ση(hydroParam, iZ, σ)
				return ση = (σ - hydroParam.σ_Min[iZ]) / (hydroParam.σ_Max[iZ] - hydroParam.σ_Min[iZ])
			end  # function: σ_2_ση
			# ------------------------------------------------------------------

end  # module: module θψ2Ks
# ............................................................