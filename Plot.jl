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

		kPa_2_mm = 101.97162 # conversion factor from kPa to mm H₂O !!!!!!!!!!!!!!!!
		
		Ψ_θΨ_Min = 0.0
		Ψ_θΨ_Max = 3000.0
		θ_θΨ_Min = 0.0
		θ_θΨ_Max = 0.58

		for iSoil = 1:2 #N_SoilSelect

			Path = path.Plots_θΨK * "Theta_h_" *string(Id_Select[iSoil]) * ".svg"
			println(Path)
			println(Id_Select[iSoil], " ", θ_θΨ[iSoil,1:N_θΨ[iSoil]], " ", Ψ_θΨ[iSoil,1:N_θΨ[iSoil]], " ", N_θΨ[iSoil]) 
			println(Id_Select[iSoil], " ", K_KΨ[iSoil,1:N_KΨ[iSoil]], " ", Ψ_KΨ[iSoil,1:N_KΨ[iSoil]], " ", N_KΨ[iSoil]) 
			println(Id_Select[iSoil], " ", N_SoilSelect)
			println(hydro.θs[iSoil], " ", hydro.θr[iSoil], " ", hydro.Ks[iSoil], " ", hydro.σ[iSoil], " ", hydro.Ψm[iSoil], " ", hydro.θsMat[iSoil], " ", hydro.σMac[iSoil], " ", hydro.ΨmMac[iSoil])
			
			println("Plotting soil $iSoil")
			Plot_CharacUnsat = GroupPlot(2, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

			push!(Plot_CharacUnsat, Axis([
				Plots.Scatter(Ψ_θΨ[iSoil,1:N_θΨ[iSoil]]*kPa_2_mm/10.0, θ_θΨ[iSoil,1:N_θΨ[iSoil]], mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"),
				], 
				style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$\theta \ [cm^3 \ cm^{-3}]$", xmin=Ψ_θΨ_Min*kPa_2_mm/10.0, xmax=Ψ_θΨ_Max*kPa_2_mm/10.0, ymin=θ_θΨ_Min, ymax=θ_θΨ_Max, xmode="log", legendPos="south west")
			)

			push!(Plot_CharacUnsat, Axis([
				Plots.Scatter(Ψ_KΨ[iSoil,1:N_KΨ[iSoil]]*kPa_2_mm/10.0, K_KΨ[iSoil,1:N_KΨ[iSoil]] * 360.0, mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"), 
				], 
				style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$K(\psi) \ [cm \ h^{-1}]$", xmode="log", legendPos="south west")
			
			)

			save(Path, Plot_CharacUnsat)
			
		end # for iSoil
		
		return
	end  # function: HYDROPARAM

	


end  # module plot
# ............................................................