# =============================================================
#		MODULE: plot
# =============================================================
module plot

	using ...wrc, ...kunsat, ..path

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : name
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYDROPARAM(θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, N_SoilSelect, hydro)

		for iSoil = 1: N_SoilSelect
			Se = wrc.θ_2_Se(0.5, iSoil, hydro)
			println(Se)
		end
		
		return
	end  # function: name

end  # module: plot
# ............................................................