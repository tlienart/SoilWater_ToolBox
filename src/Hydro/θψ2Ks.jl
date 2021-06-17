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
		function θΨ_2_KS(hydroParam, iZ::Int64, param; IsTopsoil=1, RockFragment=0.0, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ], RockFragment_Treshold=0.4, Rtol=10^-8.0, Se_Max=1.0 )

			# Determening if we have bimodal or unimodal function
				if θs / (1.0 - RockFragment) - θsMacMat / (1.0 - RockFragment) ≥ param.hydro.θs_θsMacMat # Bimodal
					FlagBimodal = true
				else
					FlagBimodal = false
				end

			# Determine when Ks increases for increasing RockFragment	
			if RockFragment > RockFragment_Treshold
				RockFragment2 = max(2.0 * RockFragment_Treshold - RockFragment, 0.0)

				θs = (θs / (1.0 - RockFragment)) * (1.0 - RockFragment2)
				
				θsMacMat = (θsMacMat / (1.0 - RockFragment)) * (1.0 - RockFragment2)

				θr = (θr / (1.0 - RockFragment)) * (1.0 - RockFragment2)
			end				

			if IsTopsoil == 1 # <>=<>=<>=<>=<>=<>
				if FlagBimodal
					τ₁    = 5.007
					τ₂    = 0.969
					τ₃    = 0.787
					τ₁Mac = 4.734
					τ₂Mac = 0.511
					τ₃Mac = 0.041
					σMac  = 0.322
				else
					τ₁ = 5.859
					τ₂ = 0.967
					τ₃ = 0.530
				end
			elseif IsTopsoil == 0  # <>=<>=<>=<>=<>=<>
				if FlagBimodal
					τ₁    = 6.444
					τ₂    = 0.859
					τ₃    = 0.408
					τ₁Mac = 3.973
					τ₂Mac = 0.642
					τ₃Mac = 0.729
					σMac  = 1.272
				else # Unimodal  
					τ₁ = 6.484
					τ₂ = 0.854
					τ₃ = 0.316
				end
			end # IsTopsoil

			# Transformation
				T1 = 10.0 ^ - τ₁
				T2 = 2.0 * (1.0 - τ₂)
				T3 =  1.0 / (1.0 - τ₃)

			Kunsat_Uni(Se) =  T1 * ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

			if FlagBimodal
				T1Mac = 10.0 ^ - τ₁Mac 
				T2mac = 2.0 * (1.0 - τ₂Mac)
				T3mac = 1.0 / (1.0 - τ₃Mac)
				Kunsat_Bim(Se) = T1Mac * ((θs - θsMacMat) ^ T3mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2mac 
				
				Kunsat_Bimodal(Se) = Kunsat_Uni(Se) + Kunsat_Bim(Se) 
				kₛ_Model = cst.KunsatModel * QuadGK.quadgk(Se -> Kunsat_Bimodal(Se), 0.0, Se_Max; rtol=Rtol)[1]
			else
				kₛ_Model = cst.KunsatModel * QuadGK.quadgk(Se -> Kunsat_Uni(Se), 0.0, Se_Max; rtol=Rtol)[1]
			end   

			kₛ_Model = min(max(hydroParam.Ks_Min[iZ], kₛ_Model), hydroParam.Ks_Max[iZ])
		return kₛ_Model
		end  # function: KUNSAT_MODEL
		
end  # module: module θψ2Ks
# ............................................................
