# =============================================================
#		MODULE: plot
# =============================================================
module plot
	import ...wrc, ...kunsat, ..path, ..cst, ..param, ..option, ..psdThetar, ..psdFunc
	using Winston
	export HYDROPARAM, PLOT_∑INFILT,  PLOT_TREANSSTEADY, PLOT_θr, PLOT_IMP_model, PLOT_∑INFILT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROPARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDROPARAM(hydro, hydroPsd, Id_Select, K_KΨ, N_KΨ, N_Psd, N_SoilSelect, N_θΨ, θ_Rpart, θ_θΨ, Ψ_KΨ, Ψ_Rpart, Ψ_θΨ; N_Se = 500)

         θ_Sim      = Array{Float64}(undef, (N_Se))
         Kunsat_Sim = Array{Float64}(undef, (N_Se))
         θ_Sim_Psd  = Array{Float64}(undef, (N_Se))
			
			# INJECTING Ψ=0 to be seen on plot using logarithmic scale to include θs 
			# Ψ_θΨ_4plot    = replace(Ψ_θΨ, 0.0=>1.0)    # replacement in [mm]
			# Ψ_Rpart_4plot = replace(Ψ_Rpart, 0.0=>1.0) # replacement in [mm]

			Ψ_θΨ_4plot    = Ψ_θΨ
			Ψ_Rpart_4plot = Ψ_Rpart

			for iSoil = 1:N_SoilSelect
				Ψ_θΨ_Max = maximum(Ψ_θΨ[iSoil,1:N_θΨ[iSoil]]) * 2.0

				Ψ_θΨ_Min = 0.8 # [mm]

				Ψ_Sim = 10.0 .^ range(log(10.0 ^ -4.0), stop=log(Ψ_θΨ_Max), length=N_Se)

				θ_θΨ_Max = maximum(θ_θΨ[iSoil,1:N_θΨ[iSoil]]) + 0.05

				# Simulated 
				for iΨ = 1:N_Se
					θ_Sim[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iSoil, hydro)
					θ_Sim_Psd[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iSoil, hydroPsd)
					if option.hydro.KunsatΨ
						Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(Ψ_Sim[iΨ], iSoil, hydro)	
					end	
				end
					
				MultiPlots = Winston.Table(1,2)
				
				# Plot Ψ(θ)
				Plot_θ_Ψ = Winston.FramedPlot(aspect_ratio=1)
					Winston.setattr(Plot_θ_Ψ.x1, label="Ψ [cm]", range=(Ψ_θΨ_Min*cst.mm_2_cm, Ψ_θΨ_Max*cst.mm_2_cm), log=true)
					Winston.setattr(Plot_θ_Ψ.y1, label="θ [cm^{3} cm^{-3}]", range=(0.0, θ_θΨ_Max))
					
					Obs_θ_Ψ = Winston.Points(Ψ_θΨ_4plot[iSoil,1:N_θΨ[iSoil]] .* cst.mm_2_cm, θ_θΨ[iSoil,1:N_θΨ[iSoil]], color="red", kind="square", size=1.5)
					Winston.setattr(Obs_θ_Ψ, label="Obs")
					
					Sim_θ_Ψ = Winston.Curve(Ψ_Sim .* cst.mm_2_cm, θ_Sim, color="red", linewidth=5)
					Winston.setattr(Sim_θ_Ψ, label="Sim")

					if option.psd.Plot_Psd_θ_Ψ && option.Psd
						Psd_θ_Ψ = Winston.Points(Ψ_Rpart_4plot[iSoil,1:N_Psd[iSoil]] .* cst.mm_2_cm, θ_Rpart[iSoil,1:N_Psd[iSoil]], color="blue", kind="circle", size=1.5)
						Winston.setattr(Psd_θ_Ψ, label="Psd")
						Sim_Psd_θ_Ψ = Winston.Curve(Ψ_Sim .* cst.mm_2_cm, θ_Sim_Psd, color="blue", linewidth=5)
						Winston.setattr(Sim_Psd_θ_Ψ, label="Sim Psd")
						legend_θ_Ψ = Winston.Legend(0.1, 0.25, [Obs_θ_Ψ, Sim_θ_Ψ, Psd_θ_Ψ, Sim_Psd_θ_Ψ])
						θ_Ψ = Winston.add(Plot_θ_Ψ, Obs_θ_Ψ, Sim_θ_Ψ, Psd_θ_Ψ, Sim_Psd_θ_Ψ, legend_θ_Ψ)
					else
						legend_θ_Ψ = Winston.Legend(0.1, 0.1, [Obs_θ_Ψ, Sim_θ_Ψ])
						θ_Ψ = Winston.add(Plot_θ_Ψ, Obs_θ_Ψ, Sim_θ_Ψ, legend_θ_Ψ)
					end

					MultiPlots[1,1] = Plot_θ_Ψ
				
				# Plot K(Ψ) 
			 	if option.hydro.KunsatΨ
				 Plot_K_Ψ = Winston.FramedPlot(aspect_ratio=1)  
					Winston.setattr(Plot_K_Ψ.x1, label="Ψ [cm]", range=(Ψ_θΨ_Min*cst.mm_2_cm, Ψ_θΨ_Max*cst.mm_2_cm), log=true)
					Winston.setattr(Plot_K_Ψ.y1, label="K(Ψ) [cm h^{-1}]")
					
					Obs_K_Ψ = Winston.Points(Ψ_KΨ[iSoil,1:N_KΨ[iSoil]] .* cst.mm_2_cm, K_KΨ[iSoil,1:N_KΨ[iSoil]] * cst.mms_2_cmh, color="red", kind="square", size=1.5)
					Winston.setattr(Obs_K_Ψ, label="Obs")

					Sim_K_Ψ = Winston.Curve(Ψ_Sim .* cst.mm_2_cm, Kunsat_Sim .* cst.mms_2_cmh, color="red", linewidth=5)
					Winston.setattr(Sim_K_Ψ, label="Sim")

					legend_K_Ψ = Winston.Legend(0.8, 0.9, [Obs_K_Ψ, Sim_K_Ψ])
					K_θ = Winston.add(Plot_K_Ψ, Obs_K_Ψ, Sim_K_Ψ, legend_K_Ψ)

					MultiPlots[1,2] = Plot_K_Ψ
				end # if option.hydro.KunsatΨ

				 Path = path.Plots_θΨK * "Lab_ThetaH_" * option.hydro.HydroModel * "_" *string(Id_Select[iSoil]) * ".svg"     
				 Winston.savefig(MultiPlots, Path)
			end # for iSoil
			
			return
		end  # function: HYDROPARAM


		function PLOT_θr(∑Psd, N_SoilSelect, hydro, paramPsd)	
			# Sorting ascending order with clay fraction
			# Array = zeros(Float64, 3, length(∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size]))
			Array = zeros(Float64, (3, N_SoilSelect))
			Array[1,:] =∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size] # Clay fraction
			Array[2,:] = paramPsd.θr_Psd[1:N_SoilSelect]
			Array[3,:] = hydro.θr[1:N_SoilSelect]
			Array = sortslices(Array, dims=2)
			Clay = Array[1,:] # Clay fraction
			θr_Psd = Array[2,:]
			θr = Array[3,:]
			
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
				θr_Sim = Winston.Points(Clay, θr, color="violet")   
				Winston.setattr(θr_Sim, label="θ_{r}")
				θrPsd = Winston.Curve(Clay, θr_Psd, color="cyan")
				Winston.setattr(θrPsd, label="θ_{r psd}")
				legend_θr_Clay = Winston.Legend(0.8, 0.15, [θr_Sim, θrPsd])
				θr_Clay = Winston.add(Plot_θr_Clay, θr_Sim, θrPsd, legend_θr_Clay) 
				MultiPlots[1,1] = Plot_θr_Clay
		   
			# Plot θr_Psd(θr)
			Plot_θrPsd_θrSim = Winston.FramedPlot(aspect_ratio=1)  
				Winston.setattr(Plot_θrPsd_θrSim.x1, label="θ_{r} [cm^3 cm^{-3}]", range=(θr_Min, θr_Max))
				Winston.setattr(Plot_θrPsd_θrSim.y1, label="θ_{r psd} [cm^3 cm^{-3}]", range=(θr_Min, θr_Max))
				θrPsd_θrSim = Winston.Points(θr, θr_Psd, color="teal")
				s = Winston.Slope(1, (0,0), kind="dotted")
				θrPsd_θrSim = Winston.add(Plot_θrPsd_θrSim, θrPsd_θrSim, s)
				MultiPlots[1,2] = Plot_θrPsd_θrSim
			
			Path = path.Plots_Psd_ThetaR
			Winston.savefig(MultiPlots, Path)
			return
		end # function: PLOT_θr



		function PLOT_IMP_model(Id_Select, Rpart, N_Psd, ∑Psd, Psd, N_SoilSelect, hydro, paramPsd)					
			for iSoil = 1:N_SoilSelect
				Rpart_Min = minimum(Rpart[iSoil,1:N_Psd[iSoil]])
				Rpart_Max = maximum(Rpart[iSoil,1:N_Psd[iSoil]])
				∑Psd_Min  = minimum(∑Psd[iSoil,1:N_Psd[iSoil]])
				∑Psd_Max  = maximum(∑Psd[iSoil,1:N_Psd[iSoil]])

				
				IntergranularMixing = zeros(Float64, N_Psd[iSoil])
				
				for iRpart = 1:N_Psd[iSoil]
					ξ2 = psdFunc.imp.∑PSD_2_ξ2(∑Psd[iSoil, iRpart]; ∑Psd_2_ξ2_β1=param.psd.imp.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=param.psd.imp.∑Psd_2_ξ2_β2)

					IntergranularMixing[iRpart] = psdFunc.imp.INTERGRANULARMIXING(Rpart[iSoil,iRpart], paramPsd.ξ1[iSoil], ξ2)
				end

				MultiPlots = Winston.Table(1,3)
				
				
				Plot_∑Psd_Rpart = Winston.FramedPlot(aspect_ratio=1)
					Winston.setattr(Plot_∑Psd_Rpart.x1, label="R_{part} [mm]", log=true) #range=(Rpart_Min, Rpart_Max), log=true)
					Winston.setattr(Plot_∑Psd_Rpart.y1, label="∑PSD") #, range=(∑Psd_Min, ∑Psd_Max))
					
					∑Psd_Rpart_points = Winston.Points(Rpart[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], color="teal", kind="square", size=1.5)
					
					∑Psd_Rpart_curve = Winston.Curve(Rpart[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], color="teal", linewidth=5) #TODO try to make it smooth 
					
					∑Psd_Rpart = Winston.add(Plot_∑Psd_Rpart, ∑Psd_Rpart_points, ∑Psd_Rpart_curve)
					MultiPlots[1,1] = Plot_∑Psd_Rpart

				Plot_Psd_Rpart = Winston.FramedPlot(aspect_ratio=1)
					Winston.setattr(Plot_Psd_Rpart.x1, label="R_{part} [mm]", log=true) #range=(Rpart_Min, Rpart_Max), log=true)
					Winston.setattr(Plot_Psd_Rpart.y1, label="PSD", range=(0.0, 0.5))
					Psd_Rpart_points = Winston.Points(Rpart[iSoil,1:N_Psd[iSoil]], Psd[iSoil,1:N_Psd[iSoil]], color="cyan", kind="square", size=1.5)
					Psd_Rpart_curve = Winston.Curve(Rpart[iSoil,1:N_Psd[iSoil]], Psd[iSoil,1:N_Psd[iSoil]], color="cyan", linewidth=5) #TODO try to make it smooth 

					Psd_Rpart = Winston.add(Plot_Psd_Rpart, Psd_Rpart_points, Psd_Rpart_curve)
				
					MultiPlots[1,2] = Plot_Psd_Rpart
				
				# Plot_NormMixing_Rpart = Winston.FramedPlot(aspect_ratio=1)
				# 	Winston.setattr(Plot_NormMixing_Rpart.x1, label="R_{part} [mm]", range=(0.0005, 0.5), log=true)
				# 	Winston.setattr(Plot_NormMixing_Rpart.y1, label="Normalised R_{part}^{-ξ(R_{Part})}", range=(0.0, 0.5))
				# 	NormMixing_Rpart_line = Winston.Points(Rpart[iSoil,1:N_Psd[iSoil]], IntergranularMixing[1:N_Psd[iSoil]], color="purple", kind="-", size=1.5)

				# 	NormMixing_Rpart = Winston.add(Plot_NormMixing_Rpart, NormMixing_Rpart_line)

				# 	MultiPlots[1,3] = Plot_NormMixing_Rpart

				Plot_NormMixing_Rpart = Winston.FramedPlot(aspect_ratio=1)
					Winston.setattr(Plot_NormMixing_Rpart.x1, label="R_{part} [mm]", log=true) #range=(Rpart_Min, Rpart_Max), log=true)
					Winston.setattr(Plot_NormMixing_Rpart.y1, label="Intergranular")
					
					Psd_Rpart_curve = Winston.Points(Rpart[iSoil,1:N_Psd[iSoil]], IntergranularMixing[1:N_Psd[iSoil]], color="cyan", linewidth=5) #TODO try to make it smooth 

					Psd_Rpart = Winston.add(Plot_NormMixing_Rpart,Psd_Rpart_curve)
				
					MultiPlots[1,3] = Plot_NormMixing_Rpart

		
				Path = path.Plots_IMP_model * "IMP_" * option.hydro.HydroModel * "_" *string(Id_Select[iSoil]) * ".svg"
				Winston.savefig(MultiPlots, Path)

			end # for iSoil
			return
		end # function: PLOT_IMP_model


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT_∑INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_∑INFILT(Id_Select, N_Infilt, N_SoilSelect, ∑Infilt_Obs, Tinfilt, ∑Infilt, infiltOutput)
		
			for iSoil=1:N_SoilSelect
				
				Plot_∑infilt_Tinfilt = Winston.FramedPlot(aspect_ratio=1)                          
					Winston.setattr(Plot_∑infilt_Tinfilt.x1, label="Time [s]")
					Winston.setattr(Plot_∑infilt_Tinfilt.y1, label="∑infiltration [mm]")

					∑infilt_Obs = Winston.Points(Tinfilt[iSoil,1:N_Infilt[iSoil]], ∑Infilt_Obs[iSoil,1:N_Infilt[iSoil]], color="violet")   
					Winston.setattr(∑infilt_Obs, label = "Obs")

					∑Infilt_Sim = Winston.Curve(Tinfilt[iSoil,1:N_Infilt[iSoil]], ∑Infilt[iSoil,1:N_Infilt[iSoil]], color="cyan")
					Winston.setattr(∑Infilt_Sim, label="Sim")

					TransSteady = Winston.Points(Tinfilt[iSoil,infiltOutput.iT_TransSteady_Data[iSoil]],∑Infilt_Obs[iSoil,infiltOutput.iT_TransSteady_Data[iSoil]], color="cyan", kind="square", size=3)
					Winston.setattr(TransSteady, label="TransSteady")
					
					legend_∑infilt_Tinfilt = Winston.Legend(0.1, 0.8, [∑infilt_Obs, ∑Infilt_Sim, TransSteady])
					∑infilt_Tinfilt = Winston.add(Plot_∑infilt_Tinfilt, ∑infilt_Obs, ∑Infilt_Sim,TransSteady, legend_∑infilt_Tinfilt)
					
					Path = path.Plots_∑infilt_Tinfilt * "INFIL_" * option.infilt.Model * "_" * option.infilt.OutputDimension *  "_" * string(Id_Select[iSoil]) *  ".svg"

				Winston.savefig(∑infilt_Tinfilt, Path)

			end # for iSoil
			
			return
		end # PLOT_∑INFILT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_TREANSSTEADY
		#		Temperorary
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_TREANSSTEADY(Id_Select, N_Infilt, N_SoilSelect, ∑Infilt_Obs, Tinfilt, ∑Infilt, infiltOutput)

			for iSoil=1:N_SoilSelect
				Plot_∑infilt_Tinfilt = Winston.FramedPlot(aspect_ratio=1)                          
					Winston.setattr(Plot_∑infilt_Tinfilt.x1, label="Time [s]")
					Winston.setattr(Plot_∑infilt_Tinfilt.y1, label="∑infiltration [mm]")

					∑infilt_Obs = Winston.Points(Tinfilt[iSoil,1:N_Infilt[iSoil]], ∑Infilt_Obs[iSoil,1:N_Infilt[iSoil]], color="violet")   
					Winston.setattr(∑infilt_Obs, label = "Obs")

					TransSteady = Winston.Points(Tinfilt[iSoil,infiltOutput.iT_TransSteady_Data[iSoil]],∑Infilt_Obs[iSoil,infiltOutput.iT_TransSteady_Data[iSoil]], color="cyan", kind="square", size=4)
					Winston.setattr(TransSteady, label="TransSteady")
					
					legend_∑infilt_Tinfilt = Winston.Legend(0.1, 0.8, [∑infilt_Obs, TransSteady])
					∑infilt_Tinfilt = Winston.add(Plot_∑infilt_Tinfilt, ∑infilt_Obs, TransSteady, legend_∑infilt_Tinfilt)
					
					Path = "C:\\JOE\\Main\\MODELS\\SOIL\\SoilWaterToolbox\\OUTPUT\\Plots\\Infiltration\\" * "INFIL_" * option.infilt.Model * "_" * option.infilt.OutputDimension *  "_" * string(Id_Select[iSoil]) *  ".svg"

				Winston.savefig(∑infilt_Tinfilt, Path)
			end # for iSoil
			
			return
		end  # function: PLOT_TREANSSTEADY

end  # module plot
# ............................................................