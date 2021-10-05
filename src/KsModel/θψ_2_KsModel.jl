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
		function KSMODEL(GroupBool_Select, hydro, ipGroup, KₛModel, KₛModel⍰, ksmodelτ, NiZ::Int64, optim, optimKsmodel, option, param; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], iGroup_Opt=1)

			for iZ=1:NiZ
				if Flag_RockFragment
					RockFragment₁ = RockFragment[iZ]
				else
					RockFragment₁ = 0.0
				end #@isdefined RockFragment

				# if Flag_IsTopsoil
				# 	IsTopsoil₁ = Int64(IsTopsoil[iZ])
				# else
				# 	IsTopsoil₁ = 1	# Default value				
				# end  # if: @isdefined IsTopsoil

				# To save time
				if GroupBool_Select[iZ]
					KₛModel[iZ] = θΨ_2_KSMODEL(hydro, option, ipGroup[iZ], iZ,  KₛModel⍰, ksmodelτ; RockFragment=RockFragment₁)
				end
			end # if: hydro.Ks[iZ] > eps(10.0)

		return KₛModel
		end  # function: KS_MODEL
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KUNSAT_MODELS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θΨ_2_KSMODEL(hydroParam, option, iGroup, iZ::Int64, KₛModel⍰, ksmodelτ; RockFragment=0.0, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ], τ₁Max=ksmodelτ.τ₁Max[iGroup], τ₂=ksmodelτ.τ₂[iGroup], τ₃=ksmodelτ.τ₃[iGroup], τ₁ηMin=ksmodelτ.τ₁ηMin[iGroup], τ₅=ksmodelτ.τ₅[iGroup], τ₄=ksmodelτ.τ₄[iGroup], τ₁Mac=ksmodelτ.τ₁Mac[iGroup], τ₂Mac=ksmodelτ.τ₂Mac[iGroup], τ₃Mac=ksmodelτ.τ₃Mac[iGroup], RockFragment_Treshold=0.4, Rtol=1.0E-3, Se_Max=0.9999 )

			# Determine when Ks increases for increasing RockFragment	
				if RockFragment > RockFragment_Treshold

					RockFragment2 = max(2.0 * RockFragment_Treshold - RockFragment, 0.0)

					θs = (θs / (1.0 - RockFragment)) * (1.0 - RockFragment2)
					
					θsMacMat = (θsMacMat / (1.0 - RockFragment)) * (1.0 - RockFragment2)

					θr = (θr / (1.0 - RockFragment)) * (1.0 - RockFragment2)
				end				

				# Model traditional 1 =================================================================
				"""Pollacco, J.A.P., Webb, T., McNeill, S., Hu, W., Carrick, S., Hewitt, A., Lilburne, L., 2017. Saturated hydraulic conductivity model computed from bimodal water retention curves for a range of New Zealand soils. Hydrol. Earth Syst. Sci. 21, 2725–2737. https://doi.org/10.5194/hess-21-2725-2017"""

				if option.ksModel.KₛModel⍰ =="KsModel_Traditional" # Original <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
					# Transforming matrix
						T1 =  10.0 ^ (τ₁Max / (τ₁Max - 1.0))
						# T1 =  (10.0 ^ - (1.0 / (1.0 - τ₁ ))) / 0.1
						T2 = 2.0 * (1.0 - τ₂)
						T3 = 1.0 / (1.0 - τ₃)

					# Transforming macro
						T1Mac = 10.0 ^ (τ₁Max * τ₁Mac / (1.0 - τ₁Max * τ₁Mac))
						T2Mac = 2.0 * (1.0 - τ₂ * τ₂Mac)
						T3Mac = 1.0 / (1.0 - τ₃Mac)
				
					return KₛModel = cst.KunsatModel * QuadGK.quadgk(Se -> KSMODEL_TRADITIONAL(Se, T1, T2, T3, T1Mac, T2Mac, T3Mac, θs, θsMacMat, θr, σ, Ψm, σMac, ΨmMac), 0.0, Se_Max; rtol=Rtol)[1]

						
				elseif option.ksModel.KₛModel⍰=="KsModel_New" # Original + Tσ₁ <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
					# Computting T1 as a function of σ
					Tσ₁ = τMODEL_σ(hydroParam, iZ, τ₁Max * τ₁ηMin, τ₁Max, σ; Inverse=false, τ₄=τ₄)

					# Transformation matrix
						T1 =  10.0 ^ (Tσ₁ / (Tσ₁ - 1.0))
						T2 = 2.0 * (1.0 - τ₂)
						T3 = 1.0 / (1.0 - τ₃)

					# Transformation macro
						T1Mac = 10.0 ^ (τ₁Max * τ₁Mac / (1.0 - τ₁Max * τ₁Mac))
						T2Mac = 2.0 * (1.0 - τ₂ * τ₂Mac)
						T3Mac = 1.0 / (1.0 - τ₃Mac)
			
					return KₛModel = cst.KunsatModel * QuadGK.quadgk(Se -> KSMODEL_TRADITIONAL(Se, T1, T2, T3, T1Mac, T2Mac, T3Mac, θs, θsMacMat, θr, σ, Ψm, σMac, ΨmMac), 0.0, Se_Max; rtol=Rtol)[1]

				elseif option.ksModel.KₛModel⍰=="KsModel_New2" # Original + Tσ₁ <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
					# Computting T1 as a function of σ
					Tσ₁ = τMODEL_σ2(hydroParam, iZ, τ₁Max * τ₁ηMin, τ₁Max, σ; τ₄=τ₄, σboundWater=τ₅)

					# Transformation matrix
						T1 =  10.0 ^ (Tσ₁ / (Tσ₁ - 1.0))
						T2 = 2.0 * (1.0 - τ₂)
						T3 = 1.0 / (1.0 - τ₃)

					# Transformation macro
						T1Mac = 10.0 ^ (τ₁Max * τ₁Mac / (1.0 - τ₁Max * τ₁Mac))
						T2Mac = 2.0 * (1.0 - τ₂ * τ₂Mac)
						T3Mac = 1.0 / (1.0 - τ₃Mac)
			
					return KₛModel = cst.KunsatModel * QuadGK.quadgk(Se -> KSMODEL_TRADITIONAL(Se, T1, T2, T3, T1Mac, T2Mac, T3Mac, θs, θsMacMat, θr, σ, Ψm, σMac, ΨmMac), 0.0, Se_Max; rtol=Rtol)[1]

				else
					error("option.ksModel.KₛModel⍰ = $(option.ksModel.KₛModel⍰) is not yet implemented try <KsModel_Traditional>; <KsModel_Tσ>; <KsModel_New>; <KsModel_NewSimplified> ")
					
				end  # if: Model=="Model?"
		end  # function: θΨ_2_KSMODEL 


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSMODEL_TRADITIONAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSMODEL_TRADITIONAL(Se, T1, T2, T3, T1Mac, T2Mac, T3Mac, θs, θsMacMat, θr, σ, Ψm, σMac, ΨmMac)

			Kunsat_Matrix =  T1 * ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

			Kunsat_Macro = T1Mac * ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
		return Kunsat = Kunsat_Matrix + Kunsat_Macro
		end  # function: KS_MODEL
	# ------------------------------------------------------------------

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSMODEL_NEW
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSMODEL_NEW(Se, T1, T2, T3, T4, T1Mac, T2Mac, T3Mac, θs, θsMacMat, θr, σ, Ψm, σMac, ΨmMac)
			
			Kunsat_Matrix = T1 *( ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2) ^ T4

			Kunsat_Macro = T1Mac * ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
		return Kunsat = Kunsat_Matrix + Kunsat_Macro
		end  # function: KS_MODEL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSMODEL_NEWSIMPLIFIED
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSMODEL_NEWSIMPLIFIED(Se, T1, T2, T3, T1Mac, θs, θsMacMat, θr, σ, Ψm, σMac, ΨmMac)
			Kunsat_Matrix =  T1 *( ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) )) ^ T2

			Kunsat_Macro = T1Mac * ((θs - θsMacMat) ^ 2.0) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) 
		return Kunsat = Kunsat_Matrix + Kunsat_Macro
		end  # function: KS_MODEL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  τMODEL_σSilt
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_σSilt_σ(hydroParam, iZ, Pₘᵢₙ, Pₘₐₓ, σ; Amplitude=0.5, σSilt_η=0.538, Pσ=3, Distribution⍰="Normal", Normalise=true, Invert=false)
			ση = τMODEL_σ(hydroParam, iZ, Pₘₐₓ, Pₘᵢₙ, σ)

			τσ_Dist = τMODEL_σSilt(hydroParam, iZ, Pₘₐₓ, Pₘᵢₙ, ση; σSilt_η=σSilt_η, Pσ=Pσ, Distribution⍰="Normal", Normalise=Normalise, Invert=Invert)

			τσ = τMODEL_σ(hydroParam, iZ, Pₘₐₓ, Pₘᵢₙ, σ)

		return τ = min(τσ + Amplitude * (τσ_Dist / (σSilt_η + 1.0)) , 1.0) * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
		end
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  τMODEL_σSilt
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_σSilt(hydroParam, iZ, Pₘᵢₙ, Pₘₐₓ,  σ; σSilt_η=0.538, Pσ=3.0, Distribution⍰="Normal", Normalise=true, Invert=false)

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
	#		FUNCTION : τMODEL_σ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_σ(hydroParam, iZ,  Pₘᵢₙ, Pₘₐₓ, σ; Inverse=false, τ₄=0.5)
			ση = σ_2_ση(hydroParam, iZ, σ)
			if Inverse
				return τ = (1.0 - ση) ^ τ₄  * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
			else
				return τ = ση ^ τ₄  * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
			end	
		end
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :τMODEL_σ2
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_σ2(hydroParam, iZ,  Pₘᵢₙ, Pₘₐₓ, σ; τ₄=0.5, σboundWater = 0.5)
			ση = σ_2_ση(hydroParam, iZ, σ)

		return τ = (Pₘₐₓ - Pₘᵢₙ) * min((ση / σboundWater) ^ τ₄, 1.0) + Pₘᵢₙ
		end
	#..................................................................



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : τMODEL_θsθr
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_θsθr(hydroParam, iZ,  Pₘᵢₙ, Pₘₐₓ, θs, θr, θsMacMat)
			θη = (θs - θr)
		return τ = (1.0 - θη)  * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
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