# =============================================================
#		MODULE: sorptivity
# =============================================================
module sorptivity
	import ..option
	export SORPTIVITY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		MAIM FUNCTION : SORPTIVITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SORPTIVITY(θ_Ini, iSoil::Int, hydro)
			if option.hydro.HydroModel == "Kosugi"
				return sorptivity.kg.SORPTIVITY(θ_Ini, iSoil::Int, hydro)
			elseif option.hydro.HydroModel == "vanGenuchten"
				return sorptivity.vg.SORPTIVITY(θ_Ini, iSoil::Int, hydro)
			end
		end # function SORPTIVITY

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<> 

	# =============================================================
	#		MODULE: kg
	# =============================================================
	module kg
		import ...wrc, ...diffusivity, ..option
		using QuadGK
		export SORPTIVITY

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : SORPTIVITY
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function SORPTIVITY(θ_Ini, iSoil, hydro; Rtol=10^-6.0)

				function SORPTIVITY_FUNC(θ, θ_Ini, iSoil, hydro)
					if option.infiltration.Sorptivity == "Parlange"
						return (hydro.θs[iSoil] + θ - 2.0 * θ_Ini) * diffusivity.kg.DIFFUSIVITY(θ, iSoil, hydro)
					elseif  option.infiltration.Sorptivity == "Option 2"
						return 2.0 * ( ((hydro.θs[iSoil] -  θ_Ini)^ 0.5) * (θ - θ_Ini)^ 0.5   ) * diffusivity.kg.DIFFUSIVITY(θ, iSoil, hydro)
					elseif  option.infiltration.Sorptivity == "Option 3"
						return 2.0 * (  (θ - θ_Ini)   ) * diffusivity.kg.DIFFUSIVITY(θ, iSoil, hydro)
					elseif  option.infiltration.Sorptivity == "Option 4"
						return 2.0 * (  (θ - θ_Ini) / ((θ - θ_Ini) / (hydro.θs[iSoil] - θ_Ini)) ^ (2.0 - 4.0 * π) ) * diffusivity.kg.DIFFUSIVITY(θ, iSoil, hydro)
					elseif  option.infiltration.Sorptivity == "Option 5"
						return 2.0 * (  (hydro.θs[iSoil] - θ_Ini)  ) * diffusivity.kg.DIFFUSIVITY(θ, iSoil, hydro)
					end # option.Infiltration
				end # function: SORPTIVITY_FUNC

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