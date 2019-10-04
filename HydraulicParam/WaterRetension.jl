  #= =============== Se =============== =#
module wrc
  	using ..option
	export Ψ_2_θDual, ∂θ∂Ψ, Ψ_2_SeDual, θ_2_ΨDual, θ_2_Se, Se_2_θ

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_θDual
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function Ψ_2_θDual(Ψ, iZ, hydro)
		if Ψ > 0.0
			Ψ = 0.0
		else
			Ψ = -Ψ
		end

		if option.θΨKmodel == "Kosugi"
			return θ = wrc.kg.Ψ_2_θDual(Ψ, iZ::Int, hydro)
		elseif option.θΨKmodel == "vanGenuchten"
			return θ = wrc.vg.Ψ_2_θ(Ψ, iZ::Int, hydro)
		end
	end # function Ψ_2_θDual


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_SeDual
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function Ψ_2_SeDual(Ψ, iZ, hydro)
		if Ψ > 0.0
			Ψ = 0.0
		else
			Ψ = -Ψ
		end

		if option.θΨKmodel == "Kosugi"
			return Se = wrc.kg.Ψ_2_SeDual(Ψ, iZ::Int, hydro)
		elseif option.θΨKmodel == "vanGenuchten"
			return Se = wrc.vg.Ψ_2_Se(Ψ, iZ::Int, hydro)
		end
	end # function Ψ_2_θDual


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θ_2_ΨDual
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function  θ_2_ΨDual(θ, iZ, hydro)
		if option.θΨKmodel == "Kosugi"
			return Ψ = -wrc.kg.θ_2_ΨDual(θ, iZ, hydro)
		elseif option.θΨKmodel == "vanGenuchten"
			return Ψ = -wrc.vg.θ_2_Ψ(θ, iZ, hydro)
		end
	end  # function  θ_2_ΨDual


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  ∂θ∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function ∂θ∂Ψ(Ψ, iZ::Int, hydro)
		if Ψ > 0.0
			Ψ = 0.0
		else
			Ψ = -Ψ
		end
	
		if option.θΨKmodel == "Kosugi"
			return ∂θ∂Ψ = wrc.kg.∂θ∂Ψ(Ψ, iZ::Int, hydro)
		elseif option.θΨKmodel == "vanGenuchten"
			return ∂θ∂Ψ = wrc.vg.∂θ∂Ψ(Ψ, iZ::Int, hydro)
		end
	end # function ∂θ∂Ψ



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θ_2_Se
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function θ_2_Se(θ, iZ::Int, hydro; θs=hydro.θs[iZ], θr=hydro.θr[iZ])
		Se = (θ - θr) / (θs - θr)
		return Se = max( min(Se, 1.0), 0.0)
	end # function θ_2_Se



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_θ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function Se_2_θ(Se, iZ::Int, hydro; θs=hydro.θs[iZ], θr=hydro.θr[iZ])
		θ = Se * (θs - θr) + θr
		return θ = max( min(θ, θs), θr)
	end # function Se_2_θ


	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	
	# =============================================================
	#		MODULE KOSUGI
	# =============================================================
	module kg
		import SpecialFunctions: erfc, erfcinv
		import Optim
		import ..wrc
		export Ψ_2_θDual, ∂θ∂Ψ, Ψ_2_SeDual, θ_2_ΨDual
	
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_θDual
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_θDual(Ψ, iZ, hydro; θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψm=hydro.Ψm[iZ], σ=hydro.σ[iZ], θsMat=hydro.θsMat[iZ], ΨmMac=hydro.ΨmMac[iZ], σMac=hydro.σMac[iZ])

			θ_Mat = (θsMat - θr) * 0.5 * erfc((log(Ψ / Ψm)) / (σ * sqrt(2.0))) + θr

			θ_Mac = (θs - θsMat) * 0.5 * erfc((log(Ψ / ΨmMac)) / (σMac * sqrt(2.0)))
			
			return θ = θ_Mac + θ_Mat
		end # function Ψ_2_θDual


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_SeDual
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_SeDual(Ψ, iZ, hydro; θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψm=hydro.Ψm[iZ], σ=hydro.σ[iZ], θsMat=hydro.θsMat[iZ], ΨmMac=hydro.ΨmMac[iZ], σMac=hydro.σMac[iZ])

			# Ψ must be negative
			θ = Ψ_2_θDual(Ψ, iZ, hydro; θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψm=hydro.Ψm[iZ], σ=hydro.σ[iZ], θsMat=hydro.θsMat[iZ], ΨmMac=hydro.ΨmMac[iZ], σMac=hydro.σMac[iZ])

			return Se = wrc.θ_2_Se(θ, iZ::Int, hydro)
		end # function Ψ_2_SeDual
	

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θ_2_ΨDuall
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θ_2_ΨDual(θ, iZ, hydro; θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψm=hydro.Ψm[iZ], σ=hydro.σ[iZ], θsMat=hydro.θsMat[iZ], ΨmMac=hydro.ΨmMac[iZ], σMac=hydro.σMac[iZ])

			function Of(Ψ, iZ, hydro)
				θmod = Ψ_2_θDual(Ψ, iZ, hydro)
				return OF = (θ - θmod) ^ 4.0
			end # Of

			if θs == θsMat
				Ψ = Ψm * exp(erfcinv(2.0 * (θ - θr) / (θs - θr)) * σ * sqrt(2.0))
			else 
				Optimization = Optim.optimize(Ψ -> Of(10.0 ^ Ψ, iZ, hydro), log10(0.001), log10(100000000.0), Optim.GoldenSection())
				Ψ = 10.0 ^ Optim.minimizer(Optimization)[1]
			end

			return Ψ
	   end # Se_2_ΨDual


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂θ∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂θ∂Ψ(Ψ, iZ, hydro; θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψm=hydro.Ψm[iZ], σ=hydro.σ[iZ], θsMat=hydro.θsMat[iZ], ΨmMac=hydro.ΨmMac[iZ], σMac=hydro.σMac[iZ]) 
			∂θ∂Ψ_Mat = -(θsMat - θr) * exp(-((log(Ψ / Ψm))^2.0) / (2.0 * σ^2.0)) / (Ψ* σ * sqrt(π*2.0))

			∂θ∂Ψ_Mac = -(θs - θsMat) * exp(-((log(Ψ / ΨmMac))^2.0) / (2.0 * σMac^2.0)) / (Ψ* σMac * sqrt(π*2.0))

			return ∂θ∂Ψ_Mat + ∂θ∂Ψ_Mac
		end # function ∂θ∂Ψ
	 
	end # module kg # ...............................................
 


	# ===============================================================================================
	#		MODULE VAN GENUCHTEN
	# ===============================================================================================
	module vg
		import ..wrc
		export Ψ_2_θ, ∂θ∂Ψ
	
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_θ(Ψ, iZ, hydro, θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψvg=hydro.Ψvg[iZ], N=hydro.N[iZ], Km=1.) # van Genuchten WRC
			M = 1.0 - Km / N
			Se = (1.0 + (Ψ / Ψvg) ^ N ) ^ (-M)
			return θ = wrc.Se_2_θ(Se, iZ, hydro)
		end #  function Ψ_2_θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_Se(Ψ, iZ, hydro, θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψvg=hydro.Ψvg[iZ], N=hydro.N[iZ], Km=1.) # van Genuchten WRC
			M = 1.0 - Km / N
			return Se = (1.0 + (Ψ / Ψvg) ^ N ) ^ (-M)
		end #  function Ψ_2_θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θ_2_Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θ_2_Ψ(θ, iZ, hydro, θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψvg=hydro.Ψvg[iZ], N=hydro.N[iZ], Km=1.) # van Genuchten WRC
			θ = max(min(θ, θs), θr)

			Se = wrc.θ_2_Se(θ, iZ, hydro) 

            M = 1. - Km / N
            Ψ = Ψvg * exp(log(exp(log(Se) / -M) - 1.) / N)
            return Ψ
        end #  Se_2_Ψ

		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂θ∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂θ∂Ψ(Ψ, iZ, hydro, θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψvg=hydro.Ψvg[iZ], N=hydro.N[iZ], Km=1.) 
			M = 1.0 - Km/N # van Genuchten
            return M * (θs-θr) / ( Ψvg*(1-M)) * ((Ψ/Ψvg) .^ (N * M)) .* (1.0 + (Ψ/Ψvg) .^ N) .^ (-M-1.)
        end #  function ∂θ∂Ψ

	end # module vg # ...............................................

end # module wrc # ...............................................