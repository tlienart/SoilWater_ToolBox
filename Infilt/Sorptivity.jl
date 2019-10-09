# =============================================================
#		MODULE: sorptivity
# =============================================================
module sorptivity

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
			function SORPTIVITY(iSoil, Se_Ini, hydro)
				Se_Ini =0.0
				θ_Ini = wrc.Se_2_θ(Se_Ini, iSoil, hydro)

				function SORPTIVITY_FUNC(θ, θ_Ini, iSoil, hydro)
					return (hydro.θs[iSoil] + θ - 2.0 * θ_Ini) * diffusivity.kg.DIFFUSIVITY(θ, iSoil, hydro)
				end
				
				return ( QuadGK.quadgk(θ -> SORPTIVITY_FUNC(θ, θ_Ini, iSoil, hydro), θ_Ini, hydro.θs[iSoil] )[1]) ^ 0.5  
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
		function SORPTIVITY()
		
			return
		end  # function: SORPTIVITY
		
	end  # module vg
	# ............................................................
	
end  # module: sorptivity
# ............................................................