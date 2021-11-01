using Symbolics, SpecialFunctions

function KUNSAT_DERIV()

	@variables θs, θsMacMat, θr,  Ψm, σ, KsMat, Ψ₁, ΨmMac, σMac, Ks, SeMat, SeMac, KsMat, KsMac


	# KsMat = Ks * (θsMacMat - θr) / (θs - θr)
	
	
	# SeMat = (θsMacMat - θr) / (θs - θr)

	# SeMac =(θs - θsMacMat) / (θs - θr)

	Se_Mat(Ψ₁) =  SeMat * 0.5 * erfc((log( Ψ₁ / Ψm)) / (σ * √2))

	Se_Mac(Ψ₁)  = SeMac * 0.5 * erfc((log(Ψ₁ / ΨmMac)) / (σMac * √2))

	Se(Ψ₁) = Se_Mat(Ψ₁) + Se_Mac(Ψ₁)

	# KsMat = Ks * SeMat			
	Kunsat_Mat(Ψ₁) =  KsMat * √Se(Ψ₁) * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2)) ^ 2

	# KsMac = Ks * SeMac
	Kunsat_Mac(Ψ₁) =  KsMac * √Se(Ψ₁) * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2)) ^ 2

	Kunsat(Ψ₁) = Kunsat_Mat(Ψ₁) + Kunsat_Mac(Ψ₁)


	D = Differential(Ψ₁)(Kunsat(Ψ₁))

	println("Derivatives")

	D2 = expand_derivatives(D)

	simplify(D2)

	println(D2)

end

KUNSAT_DERIV()