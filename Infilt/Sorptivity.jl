# =============================================================
#		MODULE: sorptivity
# =============================================================
module sorptivity
	import ..wrc, ..kunsat, ..option
	import QuadGK
	export SORPTIVITY, DIFFUSIVITY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SORPTIVITY
			# https://juliamath.github.io/QuadGK.jl/latest/
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SORPTIVITY(θ_Ini, iSoil, hydro; Rtol=10^-1.0) #-4 Rtol=10^-7.0

			function SORPTIVITY_FUNC(θ, θ_Ini, iSoil, hydro)

				Diffusivity = DIFFUSIVITY(θ, iSoil, hydro)

				if option.infilt.SorptivityModel == "Parlange" # <>=<>=<>=<>=<> strong non-linear behaviur for diffusivity 
					return Sorptivity = (hydro.θs[iSoil] + θ - 2.0 * θ_Ini) * Diffusivity

				elseif  option.infilt.SorptivityModel == "Option2"  # <>=<>=<>=<>=<>
					return Sorptivity = 2.0 * ( ((hydro.θs[iSoil] -  θ_Ini)^ 0.5) * (θ - θ_Ini)^ 0.5   ) * Diffusivity

				elseif  option.infilt.SorptivityModel == "Option3"  # <>=<>=<>=<>=<>
					return Sorptivity = 2.0 * (  (θ - θ_Ini)   ) * Diffusivity

				elseif  option.infilt.SorptivityModel == "Option4"  # Worst <>=<>=<>=<>=<>
					return Sorptivity = 2.0 * (  (θ - θ_Ini) / ((θ - θ_Ini) / (hydro.θs[iSoil] - θ_Ini)) ^ (2.0 - 4.0 * π) ) * Diffusivity

				elseif  option.infilt.SorptivityModel == "Option5"  # Best <>=<>=<>=<>=<>
					return Sorptivity = 2.0 * (  (hydro.θs[iSoil] - θ_Ini)  ) * Diffusivity
				end # option.infilt
			end # function: SORPTIVITY_FUNC

			return ( QuadGK.quadgk(θ -> SORPTIVITY_FUNC(θ, θ_Ini, iSoil, hydro), θ_Ini, hydro.θs[iSoil] - eps(), rtol=Rtol )[1]) ^ 0.5  
		end  # function: SORPTIVITY_MODEL


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : DIFFUSIVITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DIFFUSIVITY(θ, iSoil, hydro; Diffusivity_Min=10^-5.0) #-5		
			Kunsat = kunsat.θ_2_KUNSAT(θ, iSoil, hydro)
			Ψ = wrc.θ_2_ΨDual(θ, iSoil, hydro)
			∂θ∂Ψ = wrc.∂θ∂Ψ(Ψ, iSoil, hydro)

			if ∂θ∂Ψ >  Diffusivity_Min
				return Diffusivity = Kunsat / ∂θ∂Ψ
			else
				return Diffusivity = 0.0
			end
		end  # function: DIFFUSIVITY

	
end  # module: sorptivity
# ............................................................