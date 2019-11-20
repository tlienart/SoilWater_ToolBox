# =============================================================
#		MODULE: plot
# =============================================================
module plot
	import ...wrc, ...kunsat, ..path, ..cst, ..param, ..option, ..psdThetar
	using Winston
	export HYDROPARAM

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROPARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDROPARAM(Id_Select, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, N_SoilSelect, N_Psd, θ_Rpart, Ψ_Rpart, hydro; N_Se = 500)

			θ_Sim = Array{Float64}(undef, (N_Se))
			Kunsat_Sim = Array{Float64}(undef, (N_Se))
			for iSoil = 1:N_SoilSelect
				
				Ψ_θΨ_Max = maximum(Ψ_θΨ[iSoil,1:N_θΨ[iSoil]]) * 2.0
				Ψ_θΨ_Min = 0.9 # [mm]

				Ψ_Sim = 10.0 .^ range(log(10.0 ^ -4.0), stop=log(Ψ_θΨ_Max), length=N_Se)

				θ_θΨ_Max = maximum(θ_θΨ[iSoil,1:N_θΨ[iSoil]]) + 0.05

				# Simulated
				for iΨ = 1:N_Se
					θ_Sim[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iSoil, hydro)
					if option.hydro.KunsatΨ
						Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(Ψ_Sim[iΨ], iSoil, hydro)	
					end			
				end

				# println("		== Plotting Lab_ThetaH_ soil $iSoil ==")
				# Plot_CharacUnsat = GroupPlot(2, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

				# # Plotting Ψ(θ) 
				# 	push!(Plot_CharacUnsat, Axis([
				# 		Plots.Scatter(1.0 .+ Ψ_θΨ[iSoil,1:N_θΨ[iSoil]] .* cst.mm_2_cm, θ_θΨ[iSoil,1:N_θΨ[iSoil]], mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"),

				# 		Plots.Linear(Ψ_Sim .* cst.mm_2_cm, θ_Sim, mark="none", style="smooth, blue, very thick", legendentry=L"$Sim$"),
				# 		], 
				# 		style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$\theta \ [cm^3 \ cm^{-3}]$", xmin=0.1, xmax=Ψ_θΨ_Max*cst.mm_2_cm, ymin=0.0, ymax=θ_θΨ_Max, xmode="log", legendPos="south west")
				# 	)

				# # Plotting K(Ψ)
				# if option.hydro.KunsatΨ == true
				# 	push!(Plot_CharacUnsat, Axis([
				# 		Plots.Scatter(1.0 .+ Ψ_KΨ[iSoil,1:N_KΨ[iSoil]] .* cst.mm_2_cm, K_KΨ[iSoil,1:N_KΨ[iSoil]] * cst.mms_2_cmh, mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"), 

				# 		Plots.Linear(Ψ_Sim .* cst.mm_2_cm, Kunsat_Sim .* cst.mms_2_cmh, mark="none", style="smooth, blue, very thick", legendentry=L"$Sim$"),
				# 		], 
				# 		style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$K(\psi) \ [cm \ h^{-1}]$", xmin=0.1, xmax=Ψ_θΨ_Max * cst.mm_2_cm, xmode="log", legendPos="north east") 
				# 	)
				# end

				 Path = path.Plots_θΨK * "Lab_ThetaH_" * option.hydro.HydroModel * "_" *string(Id_Select[iSoil]) * ".svg"	
				# PGFPlots.save(Path, Plot_CharacUnsat)
				
				 MultiPlots = Winston.Table(1,2)
				
				# θ_Ψ=Winston.semilogx(1.0 .+ Ψ_θΨ[iSoil,1:N_θΨ[iSoil]] .* cst.mm_2_cm, θ_θΨ[iSoil,1:N_θΨ[iSoil]], "g^", Ψ_Sim .* cst.mm_2_cm, θ_Sim, "b-o")

				 # Plot Ψ(θ)
				Plot_θ_Ψ = Winston.FramedPlot(aspect_ratio=1)

					Winston.setattr(Plot_θ_Ψ.x1, label="Ψ [cm]", range=(Ψ_θΨ_Min*cst.mm_2_cm, Ψ_θΨ_Max*cst.mm_2_cm), log=true)
					Winston.setattr(Plot_θ_Ψ.y1, label="θ [cm^{3} cm^{-3}]", range=(0.0, θ_θΨ_Max))
					
					Obs_θ_Ψ = Winston.Points(0.1 .+ Ψ_θΨ[iSoil,1:N_θΨ[iSoil]] .* cst.mm_2_cm, θ_θΨ[iSoil,1:N_θΨ[iSoil]], color="red")
					Winston.setattr(Obs_θ_Ψ, label="Obs")
					
					Sim_θ_Ψ = Winston.Curve(Ψ_Sim .* cst.mm_2_cm, θ_Sim, color="blue")
					Winston.setattr(Sim_θ_Ψ, label="Sim")

					if option.psd.Plot_Psd_θ_Ψ
						Psd_θ_Ψ = Winston.Points(Ψ_Rpart[iSoil,1:N_Psd[iSoil]] .* cst.mm_2_cm, θ_Rpart[iSoil,1:N_Psd[iSoil]], color="green")
						Winston.setattr(Psd_θ_Ψ, label="Psd")
						legend_θ_Ψ = Winston.Legend(0.15, 0.2, [Obs_θ_Ψ, Sim_θ_Ψ, Psd_θ_Ψ])
						θ_Ψ = Winston.add(Plot_θ_Ψ, Obs_θ_Ψ, Sim_θ_Ψ, Psd_θ_Ψ, legend_θ_Ψ)
					else
						legend_θ_Ψ = Winston.Legend(0.15, 0.2, [Obs_θ_Ψ, Sim_θ_Ψ])
						θ_Ψ = Winston.add(Plot_θ_Ψ, Obs_θ_Ψ, Sim_θ_Ψ, legend_θ_Ψ)
					end

					MultiPlots[1,1] = Plot_θ_Ψ
				
				 # Plot K(Ψ) 
			 	if option.hydro.KunsatΨ
				Plot_K_Ψ = Winston.FramedPlot(aspect_ratio=1)  
					Winston.setattr(Plot_K_Ψ.x1, label="Ψ [cm]", range=(Ψ_θΨ_Min*cst.mm_2_cm, Ψ_θΨ_Max*cst.mm_2_cm), log=true)
					Winston.setattr(Plot_K_Ψ.y1, label="K(Ψ) [cm h^{-1}]")
					
					Obs_K_Ψ = Winston.Points(0.1 .+ Ψ_KΨ[iSoil,1:N_KΨ[iSoil]] .* cst.mm_2_cm, K_KΨ[iSoil,1:N_KΨ[iSoil]] * cst.mms_2_cmh, color="red")
					Winston.setattr(Obs_K_Ψ, label="Obs")

					Sim_K_Ψ = Winston.Curve(Ψ_Sim .* cst.mm_2_cm, Kunsat_Sim .* cst.mms_2_cmh, color="blue")
					Winston.setattr(Sim_K_Ψ, label="Sim")

					legend_K_Ψ = Winston.Legend(0.8, 0.9, [Obs_K_Ψ, Sim_K_Ψ])

					K_θ = Winston.add(Plot_K_Ψ, Obs_K_Ψ, Sim_K_Ψ, legend_K_Ψ)

					MultiPlots[1,2] = Plot_K_Ψ
				 end # if option.hydro.KunsatΨ

				 Winston.savefig(MultiPlots, Path)
			end # for iSoil
			
			return
		end  # function: HYDROPARAM


		function PLOT_θr(∑Psd, N_SoilSelect, hydro, psdparam)	
			# Minimum and maximum value
			θr_Min = param.hydro.θr_Min  
			θr_Max = param.hydro.θr_Max + 0.05
			Psd_Min = minimum(∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size]) - 0.05
			Psd_Max = maximum(∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size]) + 0.05

			MultiPlots = Winston.Table(1,2)
			
			# Plot θr(Clay)
			Plot_θr_Clay = Winston.FramedPlot(aspect_ratio=0.7)                          
			Winston.setattr(Plot_θr_Clay.x1, label="Clay [g g^{-1}]", range=(Psd_Min, Psd_Max))
			Winston.setattr(Plot_θr_Clay.y1, label="θ_{r} [cm^3 cm^{-3}]", range=(θr_Min, θr_Max))
			θr_Sim = Winston.Points(∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size], hydro.θr[1:N_SoilSelect], color="violet")   
			Winston.setattr(θr_Sim, label="θ_{r}")
			θrPsd = Winston.Curve(∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size], psdparam.θr_Psd[1:N_SoilSelect], color="cyan")
			
			Winston.setattr(θrPsd, label="θ_{r psd}")
			legend_θr_Clay = Winston.Legend(0.8, 0.15, [θr_Sim, θrPsd])
			θr_Clay = Winston.add(Plot_θr_Clay, θr_Sim, θrPsd, legend_θr_Clay) 
		   
			# Plot θr_Psd(θr)
			Plot_θrPsd_θrSim = Winston.FramedPlot(aspect_ratio=1)  
			Winston.setattr(Plot_θrPsd_θrSim.x1, label="θ_{r} [cm^3 cm^{-3}]", range=(θr_Min, θr_Max))
			Winston.setattr(Plot_θrPsd_θrSim.y1, label="θ_{r psd} [cm^3 cm^{-3}]", range=(θr_Min, θr_Max))
			θrPsd_θrSim = Winston.Points(hydro.θr[1:N_SoilSelect], psdparam.θr_Psd[1:N_SoilSelect], color="teal")
			s = Winston.Slope(1, (0,0), kind="dotted")
			θrPsd_θrSim = Winston.add(Plot_θrPsd_θrSim, θrPsd_θrSim, s)
			
			MultiPlots[1,1] = Plot_θr_Clay
			MultiPlots[1,2] = Plot_θrPsd_θrSim
			Path = path.Plots_Psd_ThetaR
			println("		== Plotting θr == ")
			Winston.savefig(MultiPlots, Path)
			# println(∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size])
			# println(length(∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size]))
			return
		end # function: PLOT_θr



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST_LAB_SEINIRANGE( ∑Infilt_Best_HydroObs_SeIniRange)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function BEST_LAB_SEINIRANGE(Id_Select,  ∑Infilt_Best_HydroObs_SeIniRange, N_SoilSelect, T_Best_HydroObs_SeIniRange)

		# 	for iSoil=1:N_SoilSelect
		# 		println("		== Plotting BestLab_SeIniRange_ soil= $iSoil ==")

		# 		Plot_BestLab = GroupPlot(1, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

		# 		push!(Plot_BestLab, Axis([
		# 			Plots.Linear(T_Best_HydroObs_SeIniRange[1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs_SeIniRange[iSoil, 1, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, red, very thick", legendentry=L"$SeIni=0.0$"),
			
		# 			Plots.Linear(T_Best_HydroObs_SeIniRange[1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs_SeIniRange[iSoil, 2, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, orange, very thick", legendentry=L"$SeIni=0.2$"),

		# 			Plots.Linear(T_Best_HydroObs_SeIniRange[1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs_SeIniRange[iSoil, 3, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, teal, very thick", legendentry=L"$SeIni=0.4$"),

		# 			Plots.Linear(T_Best_HydroObs_SeIniRange[1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs_SeIniRange[iSoil, 4, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, yellow, very thick", legendentry=L"$SeIni=0.6$"),

		# 			Plots.Linear(T_Best_HydroObs_SeIniRange[1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs_SeIniRange[iSoil, 5, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, blue, very thick", legendentry=L"$SeIni=0.8$"),
		# 			], 

		# 			style="width=10.4cm, height=6.4cm", xlabel=L"$Time \ [s]$", ylabel=L"$Infiltration \ [mm]$", legendPos="north west")
		# 		)
				
		# 		Path = path.Plots_BestLab_SeIniRange * "BestLab_SeIniRange_" *string(Id_Select[iSoil]) * ".svg"
		# 		save(Path, Plot_BestLab)

		# 	end # for
		# end  # function: BEST_LAB_SEINIRANGE( ∑Infilt_Best_HydroObs_SeIniRange)


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST_LAB
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function BEST_LAB(Id_Select, N_Infilt, N_SoilSelect, ∑Infilt_Best_HydroObs, Tinfilt_Best_HydroObs, Tinfilt, ∑Infilt)
		# 	for iSoil=1:N_SoilSelect
		# 		println("		== Plotting BestLab_ soil $iSoil ==")

		# 		Plot_BestLab = GroupPlot(1, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

		# 		push!(Plot_BestLab, Axis([
		# 			Plots.Linear(Tinfilt_Best_HydroObs[iSoil,1:param.infilt.Npoint_Infilt],  ∑Infilt_Best_HydroObs[iSoil, 1:param.infilt.Npoint_Infilt],  mark="none", style="smooth, red, very thick", legendentry=L"$HydroParam$"),
			
		# 			Plots.Scatter(Tinfilt[iSoil, 1:N_Infilt[iSoil]], ∑Infilt[iSoil, 1:N_Infilt[iSoil]], style="teal, very thick", onlyMarks=true, mark="o", markSize=4, legendentry=L"$Obs$"),
		# 			], 

		# 			style="width=10.4cm, height=6.4cm", xlabel=L"$Time \ [s]$", ylabel=L"$Infiltration \ [mm]$", legendPos="north west")
		# 		)

		# 		Path = path.Plots_BestLab * "BestLab_" *string(Id_Select[iSoil]) * ".svg"
		# 		save(Path, Plot_BestLab)
		# 	end  # for iSoil=1:N_SoilSelect

		# end # function BEST_LAB


end  # module plot
# ............................................................