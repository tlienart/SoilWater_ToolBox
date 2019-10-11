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

		for iSoil = 1:2 #N_SoilSelect

			Path = path.Plots_θΨK * "Theta_h_" *string(Id_Select[iSoil]) * ".svg"
			println(Path)
			println(Id_Select[iSoil], " ", θ_θΨ[iSoil,1:N_θΨ[iSoil]], " ", Ψ_θΨ[iSoil,1:N_θΨ[iSoil]], " ", N_θΨ[iSoil]) 
			println(Id_Select[iSoil], " ", K_KΨ[iSoil,1:N_KΨ[iSoil]], " ", Ψ_KΨ[iSoil,1:N_KΨ[iSoil]], " ", N_KΨ[iSoil]) 
			println(Id_Select[iSoil], " ", N_SoilSelect)
			println(hydro.θs[iSoil], " ", hydro.θr[iSoil], " ", hydro.Ks[iSoil], " ", hydro.σ[iSoil], " ", hydro.Ψm[iSoil], " ", hydro.θsMat[iSoil], " ", hydro.σMac[iSoil], " ", hydro.ΨmMac[iSoil])

		end
		
		return
	end  # function: HYDROPARAM

	


end  # module plot
# ............................................................