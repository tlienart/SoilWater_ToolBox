# =============================================================
#		MODULE: sorptivity
# =============================================================
module sorptivity
	import ..wrc, ..kunsat, ..option
	import QuadGK
	import SpecialFunctions: erfc, erfcinv
	export SORPTIVITY, DIFFUSIVITY


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SORPTIVITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function SORPTIVITY(θ_Ini, iSoil, hydro; Rtol=10^-3.0) #-4 Rtol=10^-7.0

		Se_Ini = wrc.θ_2_Se(θ_Ini, iSoil, hydro)
		Se_Sat = 1.0

		θhalf = (hydro.θs[iSoil] - θ_Ini) / 2.0

		function SORPTIVITY_θ(θ, hydro,  θ_Ini, Se_Ini)
			Se = wrc.θ_2_Se(θ, iSoil, hydro)

			return Sorptivity_θ = (2.0 * (hydro.θs[iSoil] - θ_Ini) / FLUXCONS_Se(Se, Se_Ini, Se_Sat)) * DIFFUSIVITY(θ, iSoil, hydro)
		end # SORPTIVITY_θ

		Sorptivity_θ = QuadGK.quadgk(θ -> SORPTIVITY_θ(θ, hydro, θ_Ini, Se_Ini), θ_Ini, θhalf, rtol=Rtol)[1]) 

	end 

	
	function FLUXCONS_Se(Se, Se_Ini, Se_Sat)
	
		if option.infilt.SorptivityModel == "Parlange"
			return FluxConcent = 2.0 * (Se - Se_Ini) / (Se_Sat + Se - 2.0 * Se_Ini)
		elseif  option.infilt.SorptivityModel == "Crank"
			return FluxConcent = exp(-erfcinv((Se - Se_Ini/(Se_Sat-Se_Ini))))
		elseif option.infilt.SorptivityModel == "Philip&Knight"
			return FluxConcent = (Se - Se_Ini) / (Se_Sat - Se_Ini)
		elseif option.infilt.SorptivityModel == "Brutsaert"
			return FluxConcent = ((Se - Se_Ini) / (Se_Sat - Se_Ini)^0.5)
		end # option.infilt

	end # FLUXCONS_Se



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

	function FLUX_CONCENTATRATION_Ψ()
	end



end  # module: sorptivity
# ............................................................