# =============================================================
#		module: kunsatModel
# =============================================================
module θψ2Ks
	import ..cst
	import QuadGK
	import SpecialFunctions: erfc, erfcinv

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KUNSAT_MODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	"""Pollacco, J.A.P., Webb, T., McNeill, S., Hu, W., Carrick, S., Hewitt, A., Lilburne, L., 2017. Saturated hydraulic conductivity model computed from bimodal water retention curves for a range of New Zealand soils. Hydrol. Earth Syst. Sci. 21, 2725–2737. https://doi.org/10.5194/hess-21-2725-2017"""
		function θΨ_2_KS(hydroParam, IsTopsoil, iZ::Int64, ksmodelτ; RockFragment=0.0, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ], τ₁=ksmodelτ.τ₁[IsTopsoil], τ₂=ksmodelτ.τ₂[IsTopsoil], τ₃=ksmodelτ.τ₃[IsTopsoil], τ₁Mac=ksmodelτ.τ₁Mac[IsTopsoil], τ₂Mac=ksmodelτ.τ₂Mac[IsTopsoil], τ₃Mac=ksmodelτ.τ₃Mac[IsTopsoil], RockFragment_Treshold=0.4, Rtol=10^-8.0,Model="Model1" )

			# Determine when Ks increases for increasing RockFragment	
				if RockFragment > RockFragment_Treshold
					RockFragment2 = max(2.0 * RockFragment_Treshold - RockFragment, 0.0)

					θs = (θs / (1.0 - RockFragment)) * (1.0 - RockFragment2)
					
					θsMacMat = (θsMacMat / (1.0 - RockFragment)) * (1.0 - RockFragment2)

					θr = (θr / (1.0 - RockFragment)) * (1.0 - RockFragment2)
				end				


			# Model traditional 1 =================================================================
			if Model=="Model1" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				# Transformation matrix
					T1 = 1.0 - 10.0 ^ -τ₁
					T2 = 2.0 * (1.0 - τ₂)
					T3 = 1.0 / (1.0 - τ₃)

					Kunsat_Matrix(Se) =  T1 * ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

				# Transformation macro
					T1Mac = 1.0 - 10.0 ^ - τ₁Mac 
					T2mac = 2.0 * (1.0 - τ₂Mac)
					T3mac = 1.0 / (1.0 - τ₃Mac)

					Kunsat_Macro(Se) = T1Mac * ((θs - θsMacMat) ^ T3mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2mac 
				
				Kunsat(Se) = Kunsat_Matrix(Se) + Kunsat_Macro(Se) 
				
				return kₛ_Model = cst.KunsatModel * QuadGK.quadgk(Se -> Kunsat(Se), 0.0, 1.0; rtol=Rtol)[1]
				
			elseif Model=="Model2"  # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
			nothing
			elseif Model=="Model3"  # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
			nothing
			elseif Model=="Model4"  # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
			nothing
			end  # if: Model=="Model1"

		# 	# Model Jesus a =====================================================================================

		# 	Func_Mixing_τ₁ = τ₁_Min + ∂θ∂Ψ_η * (τ₁_Max - τ₁_Min)
		# 	Func_Mixing_τ₂ = τ₂_Min + ∂θ∂Ψ_η * (τ₂_Max - τ₂_Min)
		# 	Func_Mixing_τ₃ = τ₃_Min + ∂θ∂Ψ_η * (τ₃_Max - τ₃_Min)

		# 	# Transformation
		# 		T1 = 1.0 - 10.0 ^ -τ₁
		# 		T2 = 2.0 * (1.0 - Func_Mixing_τ₂)
		# 		T3 = 1.0 / (1.0 - τ₃)
			
		# 		Kunsat_Matrix(Se) =  T1 * ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

		# 		T1Mac = 1.0 - 10.0 ^ - τ₁Mac 
		# 		T2mac = 2.0 * (1.0 - (τ₂Mac * ∂θ∂Ψ_Mac))
		# 		T3mac = 1.0 / (1.0 - τ₃Mac)

		# 		Kunsat_Macro(Se) = T1Mac * ((θs - θsMacMat) ^ T3mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2mac 
				
		# 		Kunsat(Se) = Kunsat_Matrix(Se) + Kunsat_Macro(Se) 
				
		# 		kₛ_Model = cst.KunsatModel * QuadGK.quadgk(Se -> Kunsat(Se), 0.0, Se_Max; rtol=Rtol)[1]

		# 	# Model Jesus b

		# 	kₛ_Model = min(max(hydroParam.Ks_Min[iZ], kₛ_Model), hydroParam.Ks_Max[iZ])
		# return kₛ_Model

		end  # function: θΨ_2_KS


		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# #		FUNCTION : ∂θ∂Ψ
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function ∂θ∂Ψ_η(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])

		# 	Mode_Mat = exp(Ψm - σ ^ 2)
		# return ∂θ∂Ψ_Mat = -  exp( -((log(Ψ₁ / Ψm)) ^ 2) / (2.0 * σ ^ 2)) / (Ψ₁ * σ * √(π * 2.0)) / Mode_Mat
		# end # function ∂θ∂Ψ
		
end  # module: module θψ2Ks
# ............................................................