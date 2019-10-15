# =============================================================
#		MODULE: diffusivity
# =============================================================
module diffusivity

	# =============================================================
	#		MODULE: kg
	# =============================================================
	module kg
		using ...wrc, ...kunsat
		export DIFFUSIVITY

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : DIFFUSIVITY
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DIFFUSIVITY(θ, iSoil, hydro)		
			Kunsat = kunsat.θ_2_KUNSAT(θ, iSoil, hydro)
			
			Ψ = wrc.kg.θ_2_ΨDual(θ, iSoil, hydro)
			
			∂Ψ∂θ = wrc.kg.∂Ψ∂θ(Ψ, iSoil, hydro)
					
			return Kunsat / ∂Ψ∂θ
		end  # function: DIFFUSIVITY
		
	end  # module: kg
	# ............................................................


	# =============================================================
	#		MODULE: vg
	# =============================================================
	module vg

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : DIFFUSIVITY
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DIFFUSIVITY()
		
			return
		end  # function: DIFFUSIVITY
		
	end  # module: vg
	# ............................................................
	
end  # module: diffusivity
# ............................................................
