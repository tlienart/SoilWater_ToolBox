# =============================================================
#		MODULE: plot
# =============================================================
module plot
	import ...wrc, ...kunsat, ..path, ..cst, ..param, ..option
	using PGFPlots, Winston
	export HYDROPARAM, BEST_LAB_SEINIRANGE, BEST_LAB

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROPARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDROPARAM(Id_Select, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, N_SoilSelect, hydro; N_Se = 50)

			θ_Sim = Array{Float64}(undef, (N_Se))
			Kunsat_Sim = Array{Float64}(undef, (N_Se))
			
			for iSoil = 1:N_SoilSelect
				
				Ψ_θΨ_Max = maximum(Ψ_θΨ[iSoil,1:N_θΨ[iSoil]]) * 10.0

				Ψ_Sim = 1.0 .+ 10.0 .^ range(log(10.0 ^ -2.0), stop=log(Ψ_θΨ_Max), length=N_Se)

				θ_θΨ_Max = maximum(θ_θΨ[iSoil,1:N_θΨ[iSoil]]) + 0.05

				# Simulated
				for iΨ = 1:N_Se
					θ_Sim[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iSoil, hydro)
					Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(Ψ_Sim[iΨ], iSoil, hydro)				
				end

				println("		== Plotting Lab_ThetaH_ soil $iSoil ==")
				Plot_CharacUnsat = GroupPlot(2, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

				# Plotting Ψ(θ) 
					push!(Plot_CharacUnsat, Axis([
						Plots.Scatter(1.0 .+ Ψ_θΨ[iSoil,1:N_θΨ[iSoil]] .* cst.mm_2_cm, θ_θΨ[iSoil,1:N_θΨ[iSoil]], mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"),

						Plots.Linear(Ψ_Sim .* cst.mm_2_cm, θ_Sim, mark="none", style="smooth, blue, very thick", legendentry=L"$Sim$"),
						], 
						style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$\theta \ [cm^3 \ cm^{-3}]$", xmin=0.1, xmax=Ψ_θΨ_Max*cst.mm_2_cm, ymin=0.0, ymax=θ_θΨ_Max, xmode="log", legendPos="south west")
					)

				# Plotting K(Ψ)
				if option.KunsatΨ == true
					push!(Plot_CharacUnsat, Axis([
						Plots.Scatter(1.0 .+ Ψ_KΨ[iSoil,1:N_KΨ[iSoil]] .* cst.mm_2_cm, K_KΨ[iSoil,1:N_KΨ[iSoil]] * cst.mms_2_cmh, mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"), 

						Plots.Linear(Ψ_Sim .* cst.mm_2_cm, Kunsat_Sim .* cst.mms_2_cmh, mark="none", style="smooth, blue, very thick", legendentry=L"$Sim$"),
						], 
						style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$K(\psi) \ [cm \ h^{-1}]$", xmin=0.1, xmax=Ψ_θΨ_Max * cst.mm_2_cm, xmode="log", legendPos="north east") 
					)
				end

				Path = path.Plots_θΨK * "Lab_ThetaH_" *string(Id_Select[iSoil]) * ".svg"	
				PGFPlots.save(Path, Plot_CharacUnsat)
				
				# MultiPlots = Winston.Table(1,2)
				
				# # Plot Ψ(θ)
				# Plot_θ_Ψ = Winston.FramedPlot(
				# xlabel="Ψ [cm]",
				# ylabel=L"$\theta \ [cm^3 \ cm^{-3}]$"                           #"θ [cm^3 \ cm^{-3}]")
				# θ_Ψ = Winston.Curve(1.0 .+ Ψ_θΨ[iSoil,1:N_θΨ[iSoil]] .* cst.mm_2_cm, θ_θΨ[iSoil,1:N_θΨ[iSoil]], "g^", Ψ_Sim .* cst.mm_2_cm, θ_Sim, "b-o")
				# Winston.add(Plot_θ_Ψ, θ_Ψ)

				# # Plot K(Ψ)
				# Plot_K_θ = Winston.FramedPlot(
				# xlabel="Ψ [cm]",
				# ylabel="K(Ψ) [cm \ h^{-1}]")
				# K_θ = Winston.Curve(1.0 .+ Ψ_KΨ[iSoil,1:N_KΨ[iSoil]] .* cst.mm_2_cm, K_KΨ[iSoil,1:N_KΨ[iSoil]] * cst.mms_2_cmh, "g^", Ψ_Sim .* cst.mm_2_cm, Kunsat_Sim .* cst.mms_2_cmh, "b-o")
				# Winston.add(Plot_K_θ, K_θ)

				# MultiPlots[1,1] = Plot_θ_Ψ
				# MultiPlots[1,2] = Plot_K_θ

				# Winston.savefig(MultiPlots, Path)


				# a=Winston.plot(1.0 .+ Ψ_θΨ[iSoil,1:N_θΨ[iSoil]] .* cst.mm_2_cm, θ_θΨ[iSoil,1:N_θΨ[iSoil]], "g^", Ψ_Sim .* cst.mm_2_cm, θ_Sim, "b-o")				
				# #Winston.savefig(a, "Test.svg")
				
				# Winston.savefig(a, Path)


			end # for iSoil
			
			return
		end  # function: HYDROPARAM


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST_LAB_SEINIRANGE( ∑Infilt_Best_HydroObs_SeIniRange)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BEST_LAB_SEINIRANGE(Id_Select,  ∑Infilt_Best_HydroObs_SeIniRange, N_SoilSelect, T_Best_HydroObs_SeIniRange)

			for iSoil=1:N_SoilSelect
				println("		== Plotting BestLab_SeIniRange_ soil= $iSoil ==")

				Plot_BestLab = GroupPlot(1, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

				push!(Plot_BestLab, Axis([
					Plots.Linear(T_Best_HydroObs_SeIniRange[1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs_SeIniRange[iSoil, 1, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, red, very thick", legendentry=L"$SeIni=0.0$"),
			
					Plots.Linear(T_Best_HydroObs_SeIniRange[1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs_SeIniRange[iSoil, 2, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, orange, very thick", legendentry=L"$SeIni=0.2$"),

					Plots.Linear(T_Best_HydroObs_SeIniRange[1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs_SeIniRange[iSoil, 3, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, teal, very thick", legendentry=L"$SeIni=0.4$"),

					Plots.Linear(T_Best_HydroObs_SeIniRange[1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs_SeIniRange[iSoil, 4, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, yellow, very thick", legendentry=L"$SeIni=0.6$"),

					Plots.Linear(T_Best_HydroObs_SeIniRange[1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs_SeIniRange[iSoil, 5, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, blue, very thick", legendentry=L"$SeIni=0.8$"),
					], 

					style="width=10.4cm, height=6.4cm", xlabel=L"$Time \ [s]$", ylabel=L"$Infiltration \ [mm]$", legendPos="north west")
				)
				
				Path = path.Plots_BestLab_SeIniRange * "BestLab_SeIniRange_" *string(Id_Select[iSoil]) * ".svg"
				save(Path, Plot_BestLab)

			end # for
		end  # function: BEST_LAB_SEINIRANGE( ∑Infilt_Best_HydroObs_SeIniRange)


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST_LAB
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BEST_LAB(Id_Select, N_Infilt, N_SoilSelect, ∑Infilt_Best_HydroObs, Tinfilt_Best_HydroObs, Tinfilt, ∑Infilt)
			for iSoil=1:N_SoilSelect
				println("		== Plotting BestLab_ soil $iSoil ==")

				Plot_BestLab = GroupPlot(1, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

				push!(Plot_BestLab, Axis([
					Plots.Linear(Tinfilt_Best_HydroObs[iSoil,1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs[iSoil, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, red, very thick", legendentry=L"$HydroParam$"),
			
					Plots.Scatter(Tinfilt[iSoil, 1:N_Infilt[iSoil]], ∑Infilt[iSoil, 1:N_Infilt[iSoil]], style="teal, very thick", onlyMarks=true, mark="o", markSize=4, legendentry=L"$Obs$"),
					], 

					style="width=10.4cm, height=6.4cm", xlabel=L"$Time \ [s]$", ylabel=L"$Infiltration \ [mm]$", legendPos="north west")
				)

				Path = path.Plots_BestLab * "BestLab_" *string(Id_Select[iSoil]) * ".svg"
				save(Path, Plot_BestLab)
			end  # for iSoil=1:N_SoilSelect

		end # function BEST_LAB


end  # module plot
# ............................................................