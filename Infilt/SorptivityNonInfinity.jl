# =============================================================
#		MODULE: sorptivity
# =============================================================
module sorptivity
	import ..wrc, ..kunsat, ..option
	import QuadGK
	import SpecialFunctions: erfc, erfcinv
	export SORPTIVITY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SORPTIVITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SORPTIVITY(θ_Ini, iSoil, hydro; ε=0.001,  Rtol=10^-2.0)

			θhalf = (hydro.θs[iSoil] + θ_Ini) / 2.0
			Ψ_Sat = 0.0
			Ψhalf = wrc.θ_2_ΨDual(θhalf, iSoil, hydro)

			function FLUXCONS_θ(θ, θ_Ini, θs)
				if option.infilt.SorptivityModel == "Parlange" # <>=<>=<>=<>=<>
					FluxConcent = 2.0 * (θ - θ_Ini) / (θs + θ - 2.0 * θ_Ini)
					return FluxConcent
		
				elseif  option.infilt.SorptivityModel == "Crank" # <>=<>=<>=<>=<>
					return FluxConcent = exp(- erfcinv( (θ - θ_Ini) / (θs - θ_Ini) ) )
		
				elseif option.infilt.SorptivityModel == "Philip&Knight" # <>=<>=<>=<>=<>
					return FluxConcent = (θ - θ_Ini) / (θs - θ_Ini)
		
				elseif option.infilt.SorptivityModel == "Brutsaert" # <>=<>=<>=<>=<>
					return FluxConcent = ((θ - θ_Ini) / (θs - θ_Ini)) ^ 0.5
				end # option.infilt.SorptivityModel
			end # FLUXCONS_θ ~~~~~~~~~~~~~~~~~


			function DIFFUSIVITY(θ, iSoil, hydro)
				Kunsat = kunsat.θ_2_KUNSAT(θ, iSoil, hydro)
				Ψ = wrc.θ_2_ΨDual(θ, iSoil, hydro)
				∂θ∂Ψ = wrc.∂θ∂Ψ(Ψ, iSoil, hydro)
				return Diffusivity = Kunsat / ∂θ∂Ψ
			end  # function: DIFFUSIVITY ~~~~~~~~~~~~~~~~~


			function SORPTIVITY_θ(θ, hydro, θ_Ini, iSoil)
				return Sorptivity_θ = DIFFUSIVITY(θ, iSoil, hydro) * (2.0 * (hydro.θs[iSoil] - θ_Ini)) / FLUXCONS_θ(θ, θ_Ini, hydro.θs[iSoil])
			end # SORPTIVITY_θ ~~~~~~~~~~~~~~~~~


			function SORPTIVITY_Ψ(Ψ, hydro, θ_Ini, iSoil)
				θ = wrc.Ψ_2_θDual(Ψ, iSoil, hydro)
				return Sorptivity_Ψ =  kunsat.θ_2_KUNSAT(θ, iSoil, hydro) * (2.0 * (hydro.θs[iSoil] - θ_Ini)) / FLUXCONS_θ(θ, θ_Ini, hydro.θs[iSoil])
			end # SORPTIVITY_Ψ ~~~~~~~~~~~~~~~~~

			Sorptivity_θ = QuadGK.quadgk(θ -> SORPTIVITY_θ(θ, hydro, θ_Ini, iSoil), θ_Ini+ε, θhalf, rtol=Rtol)[1] 

			Sorptivity_Ψ = QuadGK.quadgk(Ψ -> SORPTIVITY_Ψ(Ψ, hydro, θ_Ini, iSoil), Ψ_Sat, Ψhalf, rtol=Rtol)[1]

			return Sorptivity = (Sorptivity_θ + Sorptivity_Ψ) ^ 0.5
		end # SORPTIVITY

end  # module: sorptivity
# ............................................................