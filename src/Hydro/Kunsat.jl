module kunsat
	import ..wrc
	export Ψ_2_KUNSAT, Se_2_KUNSAT, θ_2_KUNSAT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			if  optionₘ.HydroModel⍰ == "Kosugi"
				return Kunsat = kunsat.kg.Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "Vangenuchten" ||  optionₘ.HydroModel⍰ == "VangenuchtenJules"
				return Kunsat = kunsat.vg.Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "BrooksCorey"
				return Kunsat = kunsat.bc.Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "ClappHornberger"
				return Kunsat = kunsat.ch.Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for Ψ_2_KUNSAT is not yet available")
			end
		end # function Ψ_2_KUNSAT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΨSE_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function ΨSE_2_KUNSAT(optionₘ, Ψ₁, Se, iZ::Int64, hydroParam)
			if  optionₘ.HydroModel⍰ == "Kosugi"
				return Kunsat = kunsat.kg.ΨSE_2_KUNSAT(optionₘ, Ψ₁, Se, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for ΨSE_2_KUNSAT is not yet available")
			end
		end # function Ψ_2_KUNSAT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θ_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function θ_2_KUNSAT(optionₘ, θ₁, iZ::Int64, hydroParam)
			if  optionₘ.HydroModel⍰ == "Kosugi"
				Se = wrc.θ_2_Se(θ₁, iZ, hydroParam)
				return Kunsat = Se_2_KUNSAT(optionₘ, Se, iZ, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for θ_2_KUNSAT is not yet available")
			end
		end # function θ_2_KUNSAT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam)
			Se = max(min(Se, 1.0), 0.0)

			if  optionₘ.HydroModel⍰ == "Kosugi"
				return Kunsat = kunsat.kg.Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "Vangenuchten"
				return Kunsat = kunsat.vg.Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "BrooksCorey"
				return Kunsat = kunsat.bc.Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "ClappHornberger"
				return Kunsat = kunsat.ch.Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for Se_2_KUNSAT is not yet available")
			end
		end # function Se_2_KUNSAT
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function Se_2_KR(optionₘ, Se, iZ::Int64, hydroParam)
			Se = max(min(Se, 1.0), 0.0)

			if  optionₘ.HydroModel⍰ == "Kosugi"
				return Kunsat = kunsat.kg.Se_2_KR(optionₘ, Se, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for Se_2_KR is not yet available")
			end
		end # function Se_2_KUNSAT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂K∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function ∂K∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			if  optionₘ.HydroModel⍰ == "Kosugi"
				return ∂Kunsat = kunsat.kg.∂K∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "Vangenuchten"
				return ∂Kunsat = kunsat.vg.∂K∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "BrooksCorey"
				return ∂Kunsat = kunsat.bc.∂K∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "ClappHornberger"
				return ∂Kunsat = kunsat.ch.∂K∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for ∂K∂Ψ is not yet available")
			end
		end # function ∂K∂Ψ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂K∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function ∂K∂θ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			if  optionₘ.HydroModel⍰ == "Kosugi"
				return ∂Kunsat = kunsat.kg.∂K∂θ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for ∂K∂θ is not yet available")
			end
		end


	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# =============================================================
	#		MODULE KOSUGI
	# =============================================================
	module kg
		import..wrc, ...cst
		import ForwardDiff, QuadGK
		import SpecialFunctions: erfc, erfcinv
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_KUNSAT(optionₘ, Ψ₁, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				Se = wrc.Ψ_2_SeDual(optionₘ, Ψ₁, iZ, hydroParam)

				KsMat = Ks * (θsMacMat - θr) / (θs - θr)			
				Kunsat_Mat =  KsMat * √Se * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2

				KsMac = Ks * (θs - θsMacMat) / (θs - θr)
				Kunsat_Mac =  KsMac * √Se * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2
				
				return Kunsat = Kunsat_Mat + Kunsat_Mac
			end # function Ψ_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ΨSE_2_KUNSAT(optionₘ, Ψ₁, Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				KsMat = Ks * (θsMacMat - θr) / (θs - θr)			
				Kunsat_Mat =  KsMat * √Se * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2

				KsMac = Ks * (θs - θsMacMat) / (θs - θr)
				Kunsat_Mac =  KsMac * √Se * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2
				
				return Kunsat = Kunsat_Mat + Kunsat_Mac
			end # function ΨSE_2_KUNSAT

		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Se_2_KUNSAT(optionₘ, Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				Se = max( min(Se, 1.0), 0.0)

				KsMat = Ks * (θsMacMat - θr) / (θs - θr)
				Kunsat_Mat = KsMat * √Se * (0.5 * erfc( erfcinv(2.0 * Se) + σ / √2.0 )) ^ 2

				KsMac = Ks * (θs - θsMacMat) / (θs - θr)
				Kunsat_Mac = KsMac * √Se * (0.5 * erfc( erfcinv(2.0 * Se) + σMac / √2.0 )) ^ 2

				return Kunsat = Kunsat_Mat + Kunsat_Mac
			end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Se_2_Kr
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Se_2_KR(optionₘ, Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				Se = max( min(Se, 1.0), 0.0)

				KrMat = (θsMacMat - θr) / (θs - θr)
				Kr_Mat = KrMat * √Se * (0.5*erfc( erfcinv(2.0*Se) + σ / √2.0 )) ^ 2.0

				KrMac = (θs - θsMacMat) / (θs - θr)
				Kr_Mac = KrMac * √Se * (0.5* erfc( erfcinv(2.0*Se) + σMac / √2.0 ))^2.0

				return Kr = Kr_Mat + Kr_Mac
			end # function: Se_2_KR


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂K∂Ψ numerical
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])
				
				ψ =fill(0.0::Float64, 1) 

				∂K∂Ψ_Numerical(ψ::Vector) = Ψ_2_KUNSAT(optionₘ, abs(ψ[1]), iZ, hydroParam)
				
				ψ[1] = Ψ₁
				
				Func_∂K∂Ψ_Numerical = ψ -> ForwardDiff.gradient(∂K∂Ψ_Numerical, ψ)			
				∂K∂Ψ = Func_∂K∂Ψ_Numerical(ψ)[1]
				
				if isnan(∂K∂Ψ)
					∂K∂Ψ = 0.0
				end
			return ∂K∂Ψ 
			end # function: ∂K∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂K∂Ψ analitical
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])
			
				KsMat = Ks * (θsMacMat - θr) / (θs - θr)
				∂Kunsat_Mat∂Ψ = -KsMat * (1 / (2.0*Ψ₁*√π*σ^2)) * √(erfc((log(Ψ₁ / Ψm)) / (σ*√2.0))) * erfc((log(Ψ₁ / Ψm)) / (σ*√2.0) + σ/√2.0) * exp(-0.5*((log(Ψ₁ / Ψm)) / σ + σ)))

				KsMac = Ks * (θs - θsMacMat) / (θs - θr)
				∂Kunsat_Mac∂Ψ = -KsMac * (1 / (2.0*Ψ₁*√π*σMac^2)) * √(erfc((log(Ψ₁ / ΨmMac)) / (σMac*√2.0))) * erfc((log(Ψ₁ / ΨmMac)) / (σMac*√2.0) + σMac/√2.0) * exp(-0.5*((log(Ψ₁ / ΨmMac)) / σMac + σMac)))

				∂K∂Ψ = ∂Kunsat_Mat∂Ψ + ∂Kunsat_Mac∂Ψ
			return ∂K∂Ψ 
			end # function: ∂K∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂K∂θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂θ(optionₘ, Ψ₁, Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ]) 

				P1 = 1.0 / sqrt(2.0)
				P2 = 0.125

				KsMat = Ks * (θsMacMat - θr) / (θs - θr)

				∂Kunsat_Mat∂θ = KsMat * √Se * exp(-(σ / √2.0 + erfcinv( 2.0*Se )) ^ 2 + erfcinv( 2.0*Se )^2)*erfc(σ/ √2.0 + erfcinv(2.0*Se)) + P2*KsMat*erfc(σ / √2.0 + erfcinv(2.0*Se)) ^2 / √Se

				KsMac = Ks * (θs - θsMacMat) / (θs - θr)

				∂Kunsat_Mac∂θ = KsMac*sqrt(Se)*exp(-(σMac/ √2.0 + erfcinv(2.0*Se)) ^ 2 + erfcinv(2.0*Se) ^2)*erfc(σMac / √2.0 + erfcinv(2.0*Se)) + P2*KsMac*erfc(P1*σMac + erfcinv(2.0*Se)) ^ 2 / √Se

			return ∂Kunsat_Mat∂θ + ∂Kunsat_Mac∂θ
			end #  ∂K∂θ
	end # module kg

	# =============================================================
	#		MODULE VAN GENUCHTEN
	# =============================================================
	module vg
		import ..wrc
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂Ψ

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
		#		FUNCTION vg : ∂K∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂Ψ(optionₘ, Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Ks=hydroParam.Ks[iZ], Km=hydroParam.Km[iZ], L=0.5)
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
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂Ψ

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
		#		FUNCTION bc : ∂K∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂Ψ(optionₘ, Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λbc - 2.0
		
				return ∂K∂Ψ = Ks * M / Ψbc * (Ψ₁/Ψbc) ^ (M-1.0) 

			end # function ∂K∂Ψ


	end #module bc ...............................................


	# =============================================================
	#		MODULE CLAPP AND HORNBERGER
	# =============================================================
	module ch
		import ..wrc
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂Ψ

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
		#		FUNCTION ch : ∂K∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂Ψ(optionₘ, Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λch - 2.0
		
				return ∂K∂Ψ = Ks * M / Ψch * (Ψ₁/Ψch) ^ (M-1.0) 

			end # function ∂K∂Ψ

	end #module ch ...............................................


end # module kunsat 
# ...........................................................................