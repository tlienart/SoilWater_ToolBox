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
	function HYDROPARAM(Id_Select, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, N_SoilSelect, hydro)

		for iSoil = 1:N_SoilSelect

			Path = path.Plots_θΨK * "Theta_h_" *string(Id_Select[iSoil]) * ".svg"

		end
		
		return
	end  # function: HYDROPARAM

	


end  # module plot
# ............................................................