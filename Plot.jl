# =============================================================
#		MODULE: plot
# =============================================================
module plot

	using ...wrc, ...kunsat, ..path
	using PGFPlots
	export HYDROPARAM

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROPARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYDROPARAM(θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, N_SoilSelect, hydro)

		for iSoil = 1:N_SoilSelect
			Se = wrc.θ_2_Se(0.5, iSoil, hydro)
			println(Se)
			Path = path.Plots_θΨK * "Theta_h_" *string(iSoil) * ".svg"
			println(Path)
		end
		
		return
	end  # function: HYDROPARAM

	


end  # module plot
# ............................................................