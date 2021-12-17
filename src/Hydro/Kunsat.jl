module kunsat
	import ..wrc
	export Ψ_2_KUNSAT, Se_2_KUNSAT, θ_2_KUNSAT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam)

			Ψ₁ = max(Ψ₁, 0.0)

			if  optionₘ.HydroModel⍰ == "Kosugi"
				return kunsat.kg.Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "Vangenuchten" ||  optionₘ.HydroModel⍰ == "VangenuchtenJules"
				return kunsat.vg.Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "BrooksCorey"
				return kunsat.bc.Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "ClappHornberger"
				return kunsat.ch.Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for Ψ_2_KUNSAT is not yet available")
			end
		end # function Ψ_2_KUNSAT
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΨSE_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function ΨSE_2_KUNSAT(optionₘ, Ψ₁, Se, iZ::Int64, hydroParam)

			Ψ₁ = max(Ψ₁, 0.0)

			if  optionₘ.HydroModel⍰ == "Kosugi"
				return kunsat.kg.ΨSE_2_KUNSAT(optionₘ, Ψ₁, Se, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for ΨSE_2_KUNSAT is not yet available")
			end
		end # function Ψ_2_KUNSAT
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θ_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function θ_2_KUNSAT(optionₘ, θ₁, iZ::Int64, hydroParam)
			if  optionₘ.HydroModel⍰ == "Kosugi"
				Se = wrc.θ_2_Se(θ₁, iZ, hydroParam)
				return Se_2_KUNSAT(optionₘ, Se, iZ, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for θ_2_KUNSAT is not yet available")
			end
		end # function θ_2_KUNSAT
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam)
			Se = max(min(Se, 1.0), 0.0)

			if  optionₘ.HydroModel⍰ == "Kosugi"
				return kunsat.kg.Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "Vangenuchten"
				return kunsat.vg.Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "BrooksCorey"
				return kunsat.bc.Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "ClappHornberger"
				return kunsat.ch.Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for Se_2_KUNSAT is not yet available")
			end
		end # function Se_2_KUNSAT
	#-------------------------------------------------------------------
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function Se_2_KR(optionₘ, Se, iZ::Int64, hydroParam)
			Se = max(min(Se, 1.0), 0.0)

			if  optionₘ.HydroModel⍰ == "Kosugi"
				return kunsat.kg.Se_2_KR(optionₘ, Se, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for Se_2_KR is not yet available")
			end
		end # function Se_2_KUNSAT
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂K∂ΨMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function ∂K∂ΨMODEL(optionₘ, Ψ₁, iZ::Int64, hydroParam)

			Ψ₁ = max(Ψ₁, 0.0)

			if  optionₘ.HydroModel⍰ == "Kosugi"
				return kunsat.kg.∂K∂ΨMODEL(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "Vangenuchten"
				return kunsat.vg.∂K∂ΨMODEL(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "BrooksCorey"
				return kunsat.bc.∂K∂ΨMODEL(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "ClappHornberger"
				return kunsat.ch.∂K∂ΨMODEL(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for ∂K∂ΨMODEL is not yet available")
			end
		end # function ∂K∂ΨMODEL
	#-------------------------------------------------------------------


	# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# #		FUNCTION : ∂K∂θ
	# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#   function ∂K∂Se(optionₘ, Ψ₁, Se, iZ::Int64, hydroParam)
	# 		if  optionₘ.HydroModel⍰ == "Kosugi"
	# 			return kunsat.kg.∂K∂Se(optionₘ, Ψ₁, Se, iZ::Int64, hydroParam)
	# 		else
	# 			error("$( optionₘ.HydroModel⍰) model for ∂K∂θ is not yet available")
	# 		end
	# 	end # function ∂K∂θ
	# #-------------------------------------------------------------------


	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# =============================================================
	#		MODULE KOSUGI
	# =============================================================
	module kg
		import..wrc, ...cst
		import ForwardDiff, QuadGK
		import SpecialFunctions: erfc, erfcinv
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂ΨMODEL

		θsθsMacMat =  0.001

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				Se = wrc.Ψ_2_SeDual(optionₘ, Ψ₁, iZ, hydroParam)

				KsMat = Ks * (θsMacMat - θr) / (θs - θr)			
				Kunsat_Mat =  KsMat * √Se * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0

				if θs - θsMacMat > θsθsMacMat 
					KsMac = Ks * (θs - θsMacMat) / (θs - θr)
					Kunsat_Mac =  KsMac * √Se * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0
				else
					Kunsat_Mac = 0.0
				end

			return Kunsat_Mat + Kunsat_Mac
			end # function Ψ_2_KUNSAT
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ΨSE_2_KUNSAT(optionₘ, Ψ₁, Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				KsMat = Ks * (θsMacMat - θr) / (θs - θr)			
				Kunsat_Mat =  KsMat * √Se * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0

				if θs - θsMacMat > θsθsMacMat
					KsMac = Ks * (θs - θsMacMat) / (θs - θr)
					Kunsat_Mac =  KsMac * √Se * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0
				else
					Kunsat_Mac = 0.0
				end
				
			return Kunsat_Mat + Kunsat_Mac
			end # function ΨSE_2_KUNSAT
		#-------------------------------------------------------------------

		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				Se = max( min(Se, 1.0), 0.0)

				Ψ₁ = wrc.Se_2_ΨDual(optionₘ, Se, iZ, hydroParam)

				Se_Mat = (0.5 * (θsMacMat - θr) * erfc((log( Ψ₁ / Ψm)) / (σ * √2.0))) / (θs - θr)

				KsMat = Ks * (θsMacMat - θr) / (θs - θr)
				Kunsat_Mat = KsMat * √Se * (0.5 * erfc( erfcinv(2.0 * Se_Mat) + σ / √2.0 )) ^ 2.0

				if θs - θsMacMat > θsθsMacMat
					Se_Mac = (0.5 * (θs - θsMacMat) * erfc((log(Ψ₁ / ΨmMac)) / (σMac * √2.0))) / (θs - θr)
					KsMac = Ks * (θs - θsMacMat) / (θs - θr)
					Kunsat_Mac = KsMac * √Se * (0.5 * erfc( erfcinv(2.0 * Se_Mac) + σMac / √2.0 )) ^ 2.0
				else
					Kunsat_Mac = 0.0
				end

			return Kunsat_Mat + Kunsat_Mac
			end # function Se_2_KUNSAT
		#-------------------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Se_2_Kr
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Se_2_KR(optionₘ, Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

			return Se_2_KUNSAT(optionₘ, Se, iZ, hydroParam) / Ks
			end # function: Se_2_KR
		#---------------------------------------------------------------------


		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# #		FUNCTION : ∂K∂ΨMODEL numerical
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	function ∂K∂Ψ_NUMERICAL(optionₘ, Ψ₁, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])
				
		# 		ψ =fill(0.0::Float64, 1) 

		# 		∂K∂Ψ_Numerical(ψ::Vector) = Ψ_2_KUNSAT(optionₘ, abs(ψ[1]), iZ, hydroParam)
				
		# 		ψ[1] = Ψ₁
				
		# 		Func_∂K∂Ψ_Numerical = ψ -> ForwardDiff.gradient(∂K∂Ψ_Numerical, ψ)			
		# 		∂K∂Ψ = Func_∂K∂Ψ_Numerical(ψ)[1]
				
		# 		if isnan(∂K∂Ψ)
		# 			error("∂K∂Ψ = NaN")
		# 			∂K∂Ψ = 0.0
		# 		end
		# 	return ∂K∂Ψ 
		# 	end # function: ∂K∂Ψ

			
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# #		FUNCTION : ∂K∂Ψ analitical
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	function ∂K∂Ψ_Jesus3(optionₘ, Ψ₁, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

		# 		Se = wrc.Ψ_2_SeDual(optionₘ, Ψ₁, iZ, hydroParam)

		# 		KsMat = Ks * (θsMacMat - θr) / (θs - θr)
		# 		KsMac = Ks * (θs - θsMacMat)  / (θs - θr)

		# 		Ψ₁Ψm = log(Ψ₁ / Ψm)
		# 		Ψ₁ΨmMac = log(Ψ₁ / ΨmMac)

		# 		A = Ψ₁Ψm / (σ * √2.0)
		# 		B = Ψ₁ΨmMac / (σMac * √2.0)

		# 	  ∂Kunsat_Mat∂Ψ = -(KsMat / (√(2.0*π)*Ψ₁*σ)) * (erfc(A+σ/√2.0)) * √Se * (exp(-(A+σ/√2.0)^2.0)) - (KsMat / 8*(√(2.0*π)*Ψ₁)) * ((erfc(A+σ/√2.0))^2.0) / √Se * ((θsMacMat-θr)/(θs-θr) * exp(-(A)^2.0)/σ + (θs-θsMacMat)/(θs-θr) * exp(-(B)^2.0)/σMac) 
			  
		# 	  if isnan(∂Kunsat_Mat∂Ψ)
		# 			∂Kunsat_Mat∂Ψ = 0.0
		# 		end

		# 		∂Kunsat_Mac∂Ψ = -(KsMac / (√(2.0*π)*Ψ₁*σMac)) * (erfc(B+σMac/√2.0)) * √Se * (exp(-(B+σMac/√2.0)^2.0)) - (KsMat / 8*(√(2.0*π)*Ψ₁)) * ((erfc(B+σMac/√2.0))^2.0) / √Se * ((θsMacMat-θr)/(θs-θr) * exp(-(A)^2.0)/σ + (θs-θsMacMat)/(θs-θr) * exp(-(B)^2.0)/σMac)

		# 		if isnan(∂Kunsat_Mac∂Ψ)
		# 			∂Kunsat_Mac∂Ψ = 0.0
		# 		end

		# 	return ∂Kunsat_Mat∂Ψ + ∂Kunsat_Mac∂Ψ
		# 	end # function: ∂K∂ΨMODEL analitical
		# #-----------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂K∂ΨMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂ΨMODEL(optionₘ, Ψ₁, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				if Ψ₁ < eps(100.0)
					Ψ₁ += eps(10.0)
				end

				SeMat = (θsMacMat - θr) / (θs - θr)
				SeMac =(θs - θsMacMat) / (θs - θr)

				KsMat = Ks * SeMat 
				KsMac = Ks * SeMac

				P₂ = 1.0 / √(2.0) 
				Pπ = 1.0 / √(π)
				Pc = 1.0 / √(2.0*π) 
				P8 = 1.0 / 8.0
				
				F_LOG = (P₂ * log(Ψ₁ / Ψm)) / σ
				Erfc1_Mat = erfc(P₂ * σ + F_LOG)
				ΨΨm = erfc(log(Ψ₁ / Ψm) / (√2.0 *σ))
				ΨΨmσ = exp((-(log(Ψ₁ / Ψm) / (√2.0 *σ))*log(Ψ₁ / Ψm)) / (√2.0 *σ))

				F_LOG_MAC = (P₂ * log(Ψ₁ / ΨmMac)) / σMac
				Erfc1_Mac = erfc(P₂ * σMac + F_LOG_MAC)
				ΨΨmMac = erfc(log(Ψ₁ / ΨmMac) / (√2.0 *σMac))
				ΨΨmσ_Mac = exp((-(log(Ψ₁ / ΨmMac) / (√2.0 *σMac))*log(Ψ₁ / ΨmMac)) / (√2.0 *σMac))

				∂Kunsat_Mat∂Ψ = (-Pc * KsMat*((Ψ₁ / Ψm)^-1)*Erfc1_Mat*sqrt(0.5*SeMat*ΨΨm+0.5*SeMac*ΨΨmMac)*exp(((-F_LOG) - P₂*σ) * (P₂*σ + F_LOG))) / (Ψm*σ)+ P8*KsMat*(-Pπ)*((SeMac*((Ψ₁ / ΨmMac)^-1)*ΨΨmσ_Mac) / (√2.0 *ΨmMac*σMac) + (SeMat*((Ψ₁ / Ψm)^-1)*ΨΨmσ) / (√2.0 *Ψm*σ))*(Erfc1_Mat^2)*(sqrt(0.5*SeMat*ΨΨm+0.5*SeMac*ΨΨmMac)^-1)
				
			
				if θs - θsMacMat > θsθsMacMat 
					∂Kunsat_Mac∂Ψ =  (-Pc * KsMac*((Ψ₁ / ΨmMac)^-1)*Erfc1_Mac*sqrt(0.5*SeMat*ΨΨm+0.5*SeMac*ΨΨmMac)*exp((-F_LOG_MAC - P₂*σMac)*(P₂ * σMac + F_LOG_MAC))) / (ΨmMac*σMac) + P8 *(-Pπ)*KsMac*((SeMac*((Ψ₁ / ΨmMac)^-1)*ΨΨmσ_Mac) / (√2.0 *ΨmMac*σMac) + (SeMat*((Ψ₁ / Ψm)^-1)*ΨΨmσ) / (√2.0 *Ψm*σ))*(Erfc1_Mac^2)*(sqrt(0.5*SeMat*ΨΨm+0.5*SeMac*ΨΨmMac)^-1.0)
				else
					∂Kunsat_Mac∂Ψ = 0.0
				end 

				# ∂K∂Ψ₀ = ∂Kunsat_Mac∂Ψ + ∂Kunsat_Mat∂Ψ

				# if isnan(∂K∂Ψ₀)
				# 	println(∂Kunsat_Mac∂Ψ," , " , ∂Kunsat_Mat∂Ψ)
				# 	error("∂K∂ΨMODEL =  NaN")
				# end
			return ∂Kunsat_Mac∂Ψ + ∂Kunsat_Mat∂Ψ
			end
		#--------------------------------------------------------------------------------


		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# #		FUNCTION : ∂K∂SE
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	function ∂K∂SE(optionₘ, Ψ₁, Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])
				
		# 		if Se < eps(100.0)
		# 			Se += eps(10.0)
		# 		end

		# 		P1 = 1.0 / sqrt(2.0)
		# 		P2 = 0.125

		# 		KsMat = Ks * (θsMacMat - θr) / (θs - θr)

		# 		∂Kunsat_Mat∂θ = KsMat * √Se * exp(-(σ / √2.0 + erfcinv( 2.0*Se )) ^ 2 + erfcinv(2.0*Se )^2.0)*erfc(σ/ √2.0 + erfcinv(2.0*Se)) + P2*KsMat*erfc(σ / √2.0 + erfcinv(2.0*Se)) ^2 / √Se

		# 		if θs - θsMacMat > θsθsMacMat
		# 			KsMac = Ks * (θs - θsMacMat) / (θs - θr)

		# 			∂Kunsat_Mac∂θ = KsMac*sqrt(Se)*exp(-(σMac/ √2.0 + erfcinv(2.0*Se)) ^ 2 + erfcinv(2.0*Se) ^2.0)*erfc(σMac / √2.0 + erfcinv(2.0*Se)) + P2*KsMac*erfc(P1*σMac + erfcinv(2.0*Se)) ^ 2 / √Se
		# 		else
		# 			∂Kunsat_Mac∂θ = 0.0
		# 		end

		# 	return ∂Kunsat_Mat∂θ + ∂Kunsat_Mac∂θ
		# 	end #  ∂K∂θ
		# #--------------------------------------------------------------------------------
	end # module kg



	# =============================================================
	#		MODULE VAN GENUCHTEN
	# =============================================================
	module vg
		import ..wrc
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂ΨMODEL

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Ks=hydroParam.Ks[iZ], Km=hydroParam.Km[iZ], L=0.5)
				M = 1.0 - Km / N

				Se = wrc.Ψ_2_SeDual(optionₘ, Ψ₁, iZ, hydroParam)
				return Kunsat = Ks * (Se^L) * ( 1.0 - (1.0 - Se ^ (1.0 / M) ) ^ M ) ^ 2.0
			end #function Ψ_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Se_2_KUNSAT(optionₘ, Se, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Ks=hydroParam.Ks[iZ], Km=hydroParam.Km[iZ], L=0.5)
				M = 1.0 - Km / N
			return Kunsat = Ks * (Se.^L) * ( 1. - (1. - Se.^ (1.0 / M) ) .^ M ) .^ 2.0
			end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : ∂K∂ΨMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂ΨMODEL(optionₘ, Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Ks=hydroParam.Ks[iZ], Km=hydroParam.Km[iZ], L=0.5)
				M = 1.0 - Km/N
		
				∂K∂Ψ = Ks * (L * (-M * (N * (1 / Ψvg) * (Ψ₁ / Ψvg) ^ (N - 1)) * (1.0 + (Ψ₁ / Ψvg) ^ N) ^ (-M - 1)) * ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (L - 1)) * (1.0 - (1.0 - ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ M) ^ 2.0 + Ks *
				((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ L * (2.0 * -(M * -((1.0 / M) * (-M * (N * (1 / Ψvg) * (Ψ₁ / Ψvg) ^ (N - 1)) * (1.0 + (Ψ₁ / Ψvg) ^ N) ^ (-M - 1)) * ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M - 1)) * (1.0 - ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ (M - 1)) * (1.0 - (1.0 - ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ M))

			return ∂K∂Ψ
			end # function ∂K∂Ψ

	end #module vg ...............................................


	# =============================================================
	#		MODULE BROOKS AND COOREY
	# =============================================================
	module bc
		import ..wrc
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂ΨMODEL

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λbc - 2.0

				if Ψ₁ > Ψbc
					return Kunsat = Ks * (Ψ₁ / Ψbc) ^ M 
				else 
					return Kunsat = Ks
				end

			end #function Ψ_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Se_2_KUNSAT(optionₘ, Se, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λbc - 2.0
				return Kunsat = Ks * Se .^ M 
			
			end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : ∂K∂ΨMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂ΨMODEL(optionₘ, Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λbc - 2.0
		
				return ∂K∂Ψ = Ks * M / Ψbc * (Ψ₁/Ψbc) ^ (M-1.0) 

			end # function ∂K∂Ψ


	end #module bc ...............................................


	# =============================================================
	#		MODULE CLAPP AND HORNBERGER
	# =============================================================
	module ch
		import ..wrc
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂ΨMODEL

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λch - 2.0

				if Ψ₁ > Ψch
					return Kunsat = Ks * (Ψ₁ / Ψch) ^ M 
				else 
					return Kunsat = Ks
				end

			end #function Ψ_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Se_2_KUNSAT(optionₘ, Se, iZ, hydroParam, θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψbc[iZ], λch=hydroParam.λch[iZ], Ks=hydroParam.Ks[iZ])
					
				M = -3.0 * λch - 2.0
				return Kunsat = Ks * Se .^ M 
			
			end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : ∂K∂ΨMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂ΨMODEL(optionₘ, Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λch - 2.0
		
				return ∂K∂Ψ = Ks * M / Ψch * (Ψ₁/Ψch) ^ (M-1.0) 

			end # function ∂K∂ΨMODEL

	end #module ch ...............................................


end # module kunsat 
# ...........................................................................