module kunsat
	using ..option, ..wrc
	export Ψ_2_KUNSAT, Se_2_KUNSAT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function Ψ_2_KUNSAT(Ψ, iZ::Int, hydro)
		if Ψ > 0
			Ψ = 0.0
		else
			Ψ = -Ψ
		end

		if option.θΨKmodel == "Kosugi"
			return Kunsat = kunsat.kg.Ψ_2_KUNSAT(Ψ, iZ::Int, hydro)
		elseif option.θΨKmodel == "vanGenuchten"
			return Kunsat = kunsat.vg.Ψ_2_KUNSAT(Ψ, iZ::Int, hydro)
		end
	end # function Ψ_2_KUNSAT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :Se_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function Se_2_KUNSAT(Se, iZ, hydro)
		Se = max(min(Se, 1.0), 0.0)

		if option.θΨKmodel == "Kosugi"
			return Kunsat = kunsat.kg.Se_2_KUNSAT(Se, iZ::Int, hydro)
		elseif option.θΨKmodel == "vanGenuchten"
			return Kunsat = kunsat.vg.Se_2_KUNSAT(Se, iZ::Int, hydro)
		end
	end # function Se_2_KUNSAT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂K∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function ∂K∂Ψ(Ψ, iZ, hydro)
		if Ψ > 0
			Ψ = 0.0
		else
			Ψ = -Ψ
		end
		if option.θΨKmodel == "Kosugi"
			return ∂Kunsat = kunsat.kg.∂K∂Ψ(Ψ, iZ::Int, hydro)
		elseif option.θΨKmodel == "vanGenuchten"
			return ∂Kunsat = kunsat.vg.∂K∂Ψ(Ψ, iZ::Int, hydro)
		end
	end # function ∂K∂Ψ



	# =============================================================
	#		MODULE KOSUGI
	# =============================================================
	module kg
		using ..option, ..wrc
		import SpecialFunctions: erfc, erfcinv
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_KUNSAT(Ψ, iZ, hydro; θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψm=hydro.Ψm[iZ], σ=hydro.σ[iZ], θsMat=hydro.θsMat[iZ], ΨmMac=hydro.ΨmMac[iZ], σMac=hydro.σMac[iZ], Ks=hydro.Ks[iZ])

			# Ψ must be negative
			Se = wrc.Ψ_2_SeDual(-Ψ, iZ, hydro)
			
			KsMat = Ks * (θsMat - θr) / (θs - θr)
			Kunsat_Mat =  KsMat * sqrt(Se) * (0.5 * erfc(((log(Ψ/Ψm)) / σ + σ)/sqrt(2.0)))^2.0
	
			KsMac = Ks * (θs - θsMat) / (θs - θr)
			Kunsat_Mac = KsMac * sqrt(Se) * (0.5 * erfc((( log(Ψ/ΨmMac) ) / σMac + σMac )/sqrt(2.0)))^2.0

			return Kunsat = Kunsat_Mat + Kunsat_Mac
		end # function Ψ_2_KUNSAT

		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Se_2_KUNSAT(Se, iZ, hydro; θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψm=hydro.Ψm[iZ], σ=hydro.σ[iZ], θsMat=hydro.θsMat[iZ], ΨmMac=hydro.ΨmMac[iZ], σMac=hydro.σMac[iZ], Ks=hydro.Ks[iZ])

			Se = max( min(Se, 1.0), 0.0)

			KsMat = Ks * (θsMat - θr) / (θs - θr)
			Kunsat_Mat = KsMat * sqrt(Se) * (0.5*erfc( erfcinv(2.0*Se) + σ/sqrt(2.0) )) ^2.0

			Ks = Ks * (θs - θsMat) / (θs - θr)
			Kunsat_Mac = Ks * sqrt(Se) * (0.5* erfc( erfcinv(2.0*Se) + σMac/sqrt(2.0) ))^2.0

			return Kunsat = Kunsat_Mat + Kunsat_Mac
		end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂K∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂K∂Ψ(Ψ, iZ, hydro; θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψm=hydro.Ψm[iZ], σ=hydro.σ[iZ], θsMat=hydro.θsMat[iZ], ΨmMac=hydro.ΨmMac[iZ], σMac=hydro.σMac[iZ], Ks=hydro.Ks[iZ]) 

			KsMat = Ks * (θsMat - θr) / (θs - θr)
			∂Kunsat_Mat∂Ψ = KsMat * (((0.5 * (θs - θsMat) * (((((1 / ΨmMac) * (1 / (Ψ / ΨmMac))) * (sqrt(2.0) * σMac)) / (sqrt(2.0) * σMac) ^ 2) * ((-2 * exp(-log(Ψ / ΨmMac) / (sqrt(2.0) * σMac) * (log(Ψ / ΨmMac) / (sqrt(2.0) * σMac)))) / sqrt(π))) + 0.5 * (θsMat - θr) * (((((1 / Ψm) * (1 / (Ψ / Ψm))) * (sqrt(2.0) * σ)) / (sqrt(2.0) * σ) ^ 2) * ((-2 * exp(-log(Ψ / Ψm) / (sqrt(2.0) * σ) * (log(Ψ / Ψm) / (sqrt(2.0) * σ)))) / sqrt(π)))) / (θs - θr)) * (0.5 / sqrt((0.5 * (θs - θsMat) * erfc(log(Ψ / ΨmMac) / (sqrt(2.0) * σMac)) + 0.5 * (θsMat - θr) * erfc(log(Ψ / Ψm) / (sqrt(2.0) * σ))) / (θs - θr)))) * (0.5 * erfc((log(Ψ / Ψm) / σ + σ) / sqrt(2.0))) ^ 2.0 + KsMat * sqrt((0.5 * (θs - θsMat) * erfc(log(Ψ / ΨmMac) / (sqrt(2.0) * σMac)) + 0.5 * (θsMat - θr) * erfc(log(Ψ / Ψm) / (sqrt(2.0) * σ))) / (θs - θr)) * (2.0 * (0.5 * (((((1 / Ψm) * (1 / (Ψ / Ψm))) / σ) / sqrt(2.0)) * ((-2 * exp(-(log(Ψ / Ψm) / σ + σ) / sqrt(2.0) * ((log(Ψ / Ψm) / σ + σ) / sqrt(2.0)))) / sqrt(π)))) * (0.5 * erfc((log(Ψ / Ψm) / σ + σ) / sqrt(2.0))))

			KsMac = Ks * (θs - θsMat) / (θs - θr)
			∂Kunsat_Mac∂Ψ = KsMac * (((0.5 * (θsMat - θr) * (((((1 / Ψm) * (1 / (Ψ / Ψm))) * (sqrt(2)σ)) / (sqrt(2)σ) ^ 2) * ((-2 * exp(-(log(Ψ / Ψm) / (sqrt(2)σ)) * (log(Ψ / Ψm) / (sqrt(2)σ)))) / sqrt(π))) + 0.5 * (θs - θsMat) * (((((1 / ΨmMac) * (1 / (Ψ / ΨmMac))) * (sqrt(2)σMac)) / (sqrt(2)σMac) ^ 2) * ((-2 * exp(-(log(Ψ / ΨmMac) / (sqrt(2)σMac)) * (log(Ψ / ΨmMac) / (sqrt(2)σMac)))) / sqrt(π)))) / (θs - θr)) * (0.5 / sqrt((0.5 * (θsMat - θr) * erfc(log(Ψ / Ψm) / (sqrt(2)σ)) + 0.5 * (θs - θsMat) * erfc(log(Ψ / ΨmMac) / (sqrt(2)σMac))) / (θs - θr)))) * (0.5 * erfc((log(Ψ / ΨmMac) / σMac + σMac) / sqrt(2))) ^ 2.0 + KsMac * sqrt((0.5 * (θsMat - θr) * erfc(log(Ψ / Ψm) / (sqrt(2)σ)) + 0.5 * (θs - θsMat) * erfc(log(Ψ / ΨmMac) / (sqrt(2)σMac))) / (θs - θr)) * (2.0 * (0.5 * (((((1 / ΨmMac) * (1 / (Ψ / ΨmMac))) / σMac) / sqrt(2)) * ((-2 * exp(-((log(Ψ / ΨmMac) / σMac + σMac) / sqrt(2)) * ((log(Ψ / ΨmMac) / σMac + σMac) / sqrt(2)))) / sqrt(π)))) * (0.5 * erfc((log(Ψ / ΨmMac) / σMac + σMac) / sqrt(2))))

			return ∂Kunsat = ∂Kunsat_Mac∂Ψ + ∂Kunsat_Mat∂Ψ
		end # function ∂K∂Ψ

	end # module kg 
	# =============================================================



	# =============================================================
	#		MODULE VAN GENUCHTEN
	# =============================================================
	module vg
		using ..option, ..wrc
		export Ψ_2_KUNSAT, Se_2_KUNSAT,  ∂K∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function  Ψ_2_KUNSAT(Ψ, iZ, hydro; θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψvg=hydro.Ψvg[iZ], N=hydro.N[iZ], Ks=hydro.Ks[iZ], Km=1.0, L=0.5)
			M = 1.0 - Km / N

			# Ψ must be negative
			Se = wrc.Ψ_2_SeDual(-Ψ, iZ, hydro)
		
			return Kunsat = Ks * (Se^L) * ( 1.0 - (1.0 - Se ^ (1.0 / M) ) ^ M ) ^ 2.0
		end #function Ψ_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function  Se_2_KUNSAT(Se, iZ, hydro, θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψvg=hydro.Ψvg[iZ], N=hydro.N[iZ], Ks=hydro.Ks[iZ], Km=1.0, L=0.5)
			M = 1.0 - Km/N
			return Kunsat = Ks * (Se.^L) * ( 1. - (1. - Se.^(1.0/M) ) .^ M ) .^ 2.0
		end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂K∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂K∂Ψ(Ψ, iZ, hydro, θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψvg=hydro.Ψvg[iZ], N=hydro.N[iZ], Ks=hydro.Ks[iZ], Km=1.0, L=0.5)
			M = 1.0 - Km/N
	
			∂K∂Ψ = Ks * (L * (-M * (N * (1 / Ψvg) * (Ψ / Ψvg) ^ (N - 1)) * (1.0 + (Ψ / Ψvg) ^ N) ^ (-M - 1)) * ((1.0 + (Ψ / Ψvg) ^ N) ^ -M) ^ (L - 1)) * (1.0 - (1.0 - ((1.0 + (Ψ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ M) ^ 2.0 + Ks *
			((1.0 + (Ψ / Ψvg) ^ N) ^ -M) ^ L * (2.0 * -(M * -((1.0 / M) * (-M * (N * (1 / Ψvg) * (Ψ / Ψvg) ^ (N - 1)) * (1.0 + (Ψ / Ψvg) ^ N) ^ (-M - 1)) * ((1.0 + (Ψ / Ψvg) ^ N) ^ -M) ^ (1.0 / M - 1)) * (1.0 - ((1.0 + (Ψ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ (M - 1)) * (1.0 - (1.0 - ((1.0 + (Ψ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ M))

			return ∂K∂Ψ
		end # function ∂K∂Ψ

	end #module vg ...............................................

end # module kunsat 
# ...........................................................................