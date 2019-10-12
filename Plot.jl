# =============================================================
#		MODULE: plot
# =============================================================
module plot

	using ...wrc, ...kunsat, ..path, ..cst, ..param
	using PGFPlots
	export HYDROPARAM

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROPARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYDROPARAM(Id_Select, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, N_SoilSelect, hydro; N_Se = 50)

		θ_Sim = Array{Float64}(undef, (N_Se))
		Kunsat_Sim = Array{Float64}(undef, (N_Se))
		
		for iSoil = 1:N_SoilSelect
			Path = path.Plots_θΨK * "Lab_ThetaH_" *string(Id_Select[iSoil]) * ".svg"
	
			Ψ_θΨ_Max = maximum(Ψ_θΨ[iSoil,1:N_θΨ[iSoil]]) * 10.0

			Ψ_Sim = 1.0 .+ 10.0 .^ range(log(10.0 ^ -2.0), stop=log(Ψ_θΨ_Max), length=N_Se)

			θ_θΨ_Max = maximum(θ_θΨ[iSoil,1:N_θΨ[iSoil]]) + 0.05

			# Simulated
			@simd for iΨ = 1:N_Se
				θ_Sim[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iSoil, hydro)

				Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(Ψ_Sim[iΨ], iSoil, hydro)				
			end

			println("		== Plotting soil $iSoil ==")
			Plot_CharacUnsat = GroupPlot(2, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

			# Plotting Ψ(θ) 
				push!(Plot_CharacUnsat, Axis([
					Plots.Scatter(1.0 .+ Ψ_θΨ[iSoil,1:N_θΨ[iSoil]] .* cst.mm_2_cm, θ_θΨ[iSoil,1:N_θΨ[iSoil]], mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"),

					Plots.Linear(Ψ_Sim .* cst.mm_2_cm, θ_Sim, mark="none", style="smooth, blue, very thick", legendentry=L"$Sim$"),
					], 
					style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$\theta \ [cm^3 \ cm^{-3}]$", xmin=0.1, xmax=Ψ_θΨ_Max*cst.mm_2_cm, ymin=0.0, ymax=θ_θΨ_Max, xmode="log", legendPos="south west")
				)

			# Plotting K(Ψ)
				push!(Plot_CharacUnsat, Axis([
					Plots.Scatter(1.0 .+ Ψ_KΨ[iSoil,1:N_KΨ[iSoil]] .* cst.mm_2_cm, K_KΨ[iSoil,1:N_KΨ[iSoil]] * cst.mms_2_cmh, mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"), 

					Plots.Linear(Ψ_Sim .* cst.mm_2_cm, Kunsat_Sim .* cst.mms_2_cmh, mark="none", style="smooth, blue, very thick", legendentry=L"$Sim$"),
					], 
					style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$K(\psi) \ [cm \ h^{-1}]$", xmin=0.1, xmax=Ψ_θΨ_Max * cst.mm_2_cm, xmode="log", legendPos="north east") 
				)

			save(Path, Plot_CharacUnsat)
			
		end # for iSoil
		
		return
	end  # function: HYDROPARAM

	


end  # module plot
# ............................................................