# =============================================================
#		MODULE: sorptivity
# =============================================================
module sorptivity
	using ..option
	export SORPTIVITY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		MAIM FUNCTION : SORPTIVITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SORPTIVITY(θ_Ini, iSoil::Int, hydro)
			if option.HydroModel == "Kosugi"
				return sorptivity.kg.SORPTIVITY(θ_Ini, iSoil::Int, hydro)
			elseif option.HydroModel == "vanGenuchten"
				return sorptivity.vg.SORPTIVITY(θ_Ini, iSoil::Int, hydro)
			end
		end # function SORPTIVITY

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<> 

	# =============================================================
	#		MODULE: kg
	# =============================================================
	module kg
		using ...wrc, ...diffusivity
		using QuadGK
		export SORPTIVITY

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : SORPTIVITY
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function SORPTIVITY(θ_Ini, iSoil, hydro; Rtol=10^-6.0)

				function SORPTIVITY_FUNC(θ, θ_Ini, iSoil, hydro)
					return (hydro.θs[iSoil] + θ - 2.0 * θ_Ini) * diffusivity.kg.DIFFUSIVITY(θ, iSoil, hydro)
					#return 2.0 * ( ((hydro.θs[iSoil] -  θ_Ini)^ 0.5) * (θ - θ_Ini)^ 0.5   ) * diffusivity.kg.DIFFUSIVITY(θ, iSoil, hydro)
					#return 2.0 * (  (θ - θ_Ini)   ) * diffusivity.kg.DIFFUSIVITY(θ, iSoil, hydro)
					#return 2.0 * (  (θ - θ_Ini) / ((θ - θ_Ini)/(hydro.θs[iSoil] - θ_Ini))^(2-4/pi) ) * diffusivity.kg.DIFFUSIVITY(θ, iSoil, hydro)
					#return 2.0 * (  (hydro.θs[iSoil] - θ_Ini)  ) * diffusivity.kg.DIFFUSIVITY(θ, iSoil, hydro)
				end

				return ( QuadGK.quadgk(θ -> SORPTIVITY_FUNC(θ, θ_Ini, iSoil, hydro), θ_Ini, hydro.θs[iSoil] - eps() )[1]) ^ 0.5  
			end  # function: SORPTIVITY


	end  # module: kg
	# ............................................................

	# =============================================================
	#		MODULE: vg
	# =============================================================
	module vg
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : SORPTIVITY
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SORPTIVITY(θ_Ini, iSoil, hydro)
		
			return
		end  # function: SORPTIVITY
		
	end  # module vg
	# ............................................................
	
end  # module: sorptivity
# ............................................................