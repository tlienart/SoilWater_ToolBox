using Calculus, SpecialFunctions

function KUNSAT_DERIV2()


   # θ_Mat = ((θsMacMat - θr) * 0.5 * erfc((log( Ψ₁ / Ψm)) / (σ * √2.0))) / (θs - θr)

   # Kunsat_Mat =  KsMat * ( ((θsMacMat - θr) * 0.5 * erfc((log( Ψ₁ / Ψm)) / (σ * (2.0)^0.5))) / (θs - θr))^0.5 * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / (2.0)^0.5)) ^ 2.0


   D = differentiate("KsMat * ( ((θsMacMat - θr) * 0.5 * erfc((log( Ψ₁ / Ψm)) / (σ * (2.0)^0.5))) / (θs - θr))^0.5 * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / (2.0)^0.5)) ^ 2.0", :Ψ₁)


   println(D)

end

KUNSAT_DERIV2()