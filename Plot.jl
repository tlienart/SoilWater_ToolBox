# =============================================================
#		MODULE: plot
	# julia> ] # go to Pkg mode
	# (v1.1) pkg> dev PackageCompiler
	# # go back REPL mode
	# julia > compile_package("Plots", "GR")
# =============================================================
module plot
	import ...wrc, ...kunsat, ..path, ..cst, ..param, ..option, ..psdThetar, ..psdFunc
	using Plots, Plots.PlotMeasures
	using LaTeXStrings
	export HYDROPARAM, PLOT_∑INFILT,  PLOT_TREANSSTEADY, PLOT_θr, PLOT_IMP_model, PLOT_∑INFILT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROPARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, Id_Select, hydro; N_Se = 500)

			θ_Sim      = Array{Float64}(undef, (N_Se))
			Kunsat_Sim = Array{Float64}(undef, (N_Se))
			θ_Sim_Psd  = Array{Float64}(undef, (N_Se))

			for iSoil = 1:N_SoilSelect
				Ψ_θΨ_Max = maximum(Ψ_θΨ[iSoil,1:N_θΨ[iSoil]]) * 2.0 # [mm]

				Ψ_θΨ_Min = 0.01 # [mm]

				Ψ_Sim = 10.0 .^ range(log(10.0 ^ -4.0), stop=log(Ψ_θΨ_Max), length=N_Se)

				θ_θΨ_Max = maximum(θ_θΨ[iSoil,1:N_θΨ[iSoil]]) + 0.1

				# Simulated 
				for iΨ = 1:N_Se
                    θ_Sim[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iSoil, hydro)
                    # θ_Sim_Psd[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iSoil, hydroPsd)

					if option.hydro.KunsatΨ
						Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(Ψ_Sim[iΨ], iSoil, hydro)	
					end	
				end

				# == Plot_θ_Ψ ==
				gr()
					# Plot_θ_Ψ: Observed
					X = Ψ_θΨ[iSoil,1:N_θΨ[iSoil]] .* cst.mm_2_cm
					Y = θ_θΨ[iSoil,1:N_θΨ[iSoil]]
					Plot_θ_Ψ = Plots.plot(X ,Y, seriestype=:scatter, label="Obs", color= :red, shape= :square, markersize=4)

					# Plot_θ_Ψ: Simulated
					X = Ψ_Sim[1:N_Se] .* cst.mm_2_cm
					Y = θ_Sim[1:N_Se]
					Plots.plot!(X, Y, seriestype=:line, label="Sim", color= :blue, lw=2)

					# Plot_θ_Ψ: Total porosity point
					X = zeros(Float64,1)
					X[1] = 0.1 .* cst.mm_2_cm
					Y = zeros(Float64,1)
					Y[1] = hydro.Φ[iSoil]
					Plots.plot!(X, Y, seriestype=:scatter, label= L"\phi", color= :green, shape= :square, markersize=4) 

					# Plot_θ_Ψ: General attributes
					xlabel!(L"\psi \ [cm]")
					ylabel!(L"\theta \ [cm^3 cm^{-3}]")
					Plots.plot!(xlims = (Ψ_θΨ_Min*cst.mm_2_cm*10.0, Ψ_θΨ_Max*cst.mm_2_cm), ylims = (0.0, hydro.Φ[iSoil]+0.1), xscale= :log10, size=(600,600))

					# == Plot_K_Ψ ==
					if option.hydro.KunsatΨ
						# Plot_K_Ψ: Obs K_Ψ
						X = Ψ_KΨ[iSoil,1:N_KΨ[iSoil]] .* cst.mm_2_cm
						Y = K_KΨ[iSoil,1:N_KΨ[iSoil]] * cst.mms_2_cmh
						Plot_K_Ψ = Plots.plot(X ,Y, seriestype=:scatter, label="Obs", color= :red, shape= :square, markersize=4)

						# Plot_K_Ψ: Sim K_Ψ
						X = Ψ_Sim .* cst.mm_2_cm
						Y = Kunsat_Sim .* cst.mms_2_cmh
						Plots.plot!(X, Y, seriestype=:line, label="Sim", color= :blue, lw=2)

						# Plot_K_Ψ: Ks
						X = zeros(Float64,1)
						X[1] = 0.1* cst.mm_2_cm
						Y = zeros(Float64,1)
						Y[1] =  hydro.Ks[iSoil] * cst.mms_2_cmh
						Plots.plot!(X, Y, seriestype=:scatter, label= "Ks", color= :green, shape= :square, markersize=4) 

						# General attributes
						xlabel!(L"\psi \ [cm]")
						ylabel!(L" K (\psi) \ [cm h^{-1}]")
						Plots.plot!(xlims = (Ψ_θΨ_Min*cst.mm_2_cm, Ψ_θΨ_Max*cst.mm_2_cm), ylims = (0.0, hydro.Ks[iSoil] * cst.mms_2_cmh), xscale= :log10, size=(600,600))

					end # if option.hydro.KunsatΨ

			

				# 	if option.psd.Plot_Psd_θ_Ψ && option.Psd
				# 		Psd_θ_Ψ = Winston.Points(Ψ_Rpart[iSoil,1:N_Psd[iSoil]] .* cst.mm_2_cm, θ_Rpart[iSoil,1:N_Psd[iSoil]], color="blue", kind="circle", size=1.5)


				# 		Psd_θ_Ψ = Winston.Points(Ψ_Rpart[iSoil,1:N_Psd[iSoil]] .* cst.mm_2_cm, θ_Rpart[iSoil,1:N_Psd[iSoil]], color="blue", kind="circle", size=1.5)

				# 		Winston.setattr(Psd_θ_Ψ, label="Psd")
				# 		Sim_Psd_θ_Ψ = Winston.Curve(Ψ_Sim .* cst.mm_2_cm, θ_Sim_Psd, color="green", linewidth=5)

				# 		Winston.setattr(Sim_Psd_θ_Ψ, label="Sim Psd")
				# 		legend_θ_Ψ = Winston.Legend(0.1, 0.25, [Obs_θ_Ψ, Sim_θ_Ψ, Psd_θ_Ψ, Sim_Psd_θ_Ψ])
				# 		θ_Ψ = Winston.add(Plot_θ_Ψ, Obs_θ_Ψ, Sim_θ_Ψ, Psd_θ_Ψ, Sim_Psd_θ_Ψ, Sim_Φ, legend_θ_Ψ)
				# 	else
				# 		legend_θ_Ψ = Winston.Legend(0.1, 0.5, [Obs_θ_Ψ, Sim_θ_Ψ, Sim_Φ])
				# 		θ_Ψ = Winston.add(Plot_θ_Ψ, Obs_θ_Ψ, Sim_θ_Ψ, legend_θ_Ψ, Sim_Φ)
				# 	end


				Path = path.Plots_θΨK * "Lab_ThetaH_" * option.hydro.HydroModel * "_" *string(Id_Select[iSoil]) * ".svg"     
				Plot = Plots.plot(Plot_θ_Ψ, Plot_K_Ψ)
				Plots.savefig(Plot, Path)

			end # for iSoil
				
			return
		end  # function: HYDROPARAM



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT_θr
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_θr(∑Psd, N_SoilSelect, hydro, paramPsd)	
			# Sorting ascending order with clay fraction
				Array = zeros(Float64, 3, length(∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size]))
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
			
			# PLOT 1 <>=<>=<>=<>=<>=<>
			# Plot θr(Clay)
				X = Clay
				Y = θr
				Plot_θr = Plots.plot(X ,Y, seriestype=:scatter, label=L"\theta _{r}", color= :violet, shape= :square, markersize=4)

			# Plot θr_psd(Clay)
				X = Clay
				Y = θr_Psd
				Plots.plot!(X ,Y, seriestype=:line, label=L"\theta _{r psd}", color= :cyan, lw=2)

			# General attributes
				xlabel!(L"\theta _{r} [cm^{3} \ cm^{-3}]")
				ylabel!(L"Clay \ [g \ g^{-1}]")                         
				Plots.plot!(xlims = (Psd_Min, Psd_Max), ylims = (θr_Min, θr_Max), size=(600,600))
	
			# PLOT 2 <>=<>=<>=<>=<>=<>
			# Plot θr_Psd(θr)
				X = θr
				Y = θr_Psd
				Plot_θr_Psd = Plots.plot(X ,Y, seriestype=:scatter, color=:violet, shape=:square, markersize=4)
				Plots.plot!(X, Y, seriestype=:line, label="Sim", color= :blue, lw=2)

				
				# 1:1 line
				X = θr_Min:θr_Max
				Y = X
				Plots.plot!(X, Y, seriestype=:line, label="1:1", color= :black, linestyle= :dot, lw=2)

			# General attributes
				xlabel!(L"\theta _{r} [cm^3 cm^{-3}]")
				ylabel!(L"\theta _{r psd} [cm^3 cm^{-3}]")
				Plots.plot!(xlims = (θr_Min, θr_Max), ylims = (θr_Min, θr_Max), size=(600,600))
	
			Path = path.Plots_Psd_θr
			Plot = Plots.plot(Plot_θr, Plot_θr_Psd)
			Plots.savefig(Plot, Path)
			return
		end # function: PLOT_θr


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT_IMP_model
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
				end # for iRpart = 1:N_Psd[iSoil]

				# << PLOT 1 >>
					# Plot_∑Psd_Rpart
						X = Rpart[iSoil,1:N_Psd[iSoil]]
						Y = ∑Psd[iSoil,1:N_Psd[iSoil]]
						Plot_∑Psd_Rpart = Plots.plot(X ,Y, seriestype=:scatter, label="Obs", color= :teal, shape= :square, markersize= 4)

					# Plot_∑Psd_Rpart: General attributes
						xlabel!(L"R_{part} [mm]")
						ylabel!(L"\sum \ PSD")
						
					Plots.plot!(xlims = (Rpart_Min, Rpart_Max), ylims = (∑Psd_Min, ∑Psd_Max), xscale= :log10, aspect_ratio= 1)


				# << PLOT 2 >>
					# Plot_Psd_Rpart
						X = Rpart[iSoil,1:N_Psd[iSoil]]
						Y = Psd[iSoil,1:N_Psd[iSoil]]
						Plot_∑Psd_Psd = Plots.plot(X ,Y, seriestype=:scatter, label="Obs", color= :teal, shape= :square, markersize= 4)

						xlabel!(L"R_{part} [mm]")
						ylabel!(L"PSD [mm]")
						Plots.plot!(xlims = (Rpart_Min, Rpart_Max), ylims= (0.0, 0.5), xscale= :log10, aspect_ratio= 1)


				# << PLOT 3 >>
					# Plot NormMixing_Rpart
						X = Rpart[iSoil,1:N_Psd[iSoil]]
						Y = IntergranularMixing[1:N_Psd[iSoil]]
						Plot_NormMixing_Rpart = ∑Psd_Psd = Plots.plot(X, Y, seriestype=:scatter, label="Obs", color= :teal, shape= :square, markersize=4, xscale= :log10)

						xlabel!(L"R_{part} [mm]")
						ylabel!(L"Intergranular")
						Plots.plot!(xlims= (Rpart_Min, Rpart_Max), xscale=:log10, aspect_ratio=1)

				Plot = Plots.plot(Plot_∑Psd_Rpart, Plot_∑Psd_Psd, Plot_NormMixing_Rpart, layout= (1,3))
				Path = path.Plots_IMP_model * "IMP_" * option.hydro.HydroModel * "_" *string(Id_Select[iSoil]) * ".svg"
				Plots.savefig(Plot, Path)

			end # for iSoil
			return	
		end # function: PLOT_IMP_model


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT_∑INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_∑INFILT(Id_Select, N_Infilt, N_SoilSelect, ∑Infilt_Obs, Tinfilt, ∑Infilt, infiltOutput)
		
			for iSoil=1:N_SoilSelect	
				# << PLOT 1 >>
					# Plot_∑infilt_Obs
						X = Tinfilt[iSoil,1:N_Infilt[iSoil]]
						Y = ∑Infilt_Obs[iSoil,1:N_Infilt[iSoil]]
						Plot_∑infilt_Obs = Plots.plot(X, Y, seriestype=:scatter, label="Obs", color= :violet, shape= :square, markersize=4) 

					# Plot_∑infilt_Sim
						X = Tinfilt[iSoil,1:N_Infilt[iSoil]]
						Y = ∑Infilt[iSoil,1:N_Infilt[iSoil]]
						Plots.plot!(X, Y, seriestype=:line, label="Sim", color= :cyan, shape= :square, markersize=4) 

					# TransSteady
						X = zeros(Float64,1)
						Y = zeros(Float64,1)
						X[1] = Tinfilt[iSoil,infiltOutput.iT_TransSteady_Data[iSoil]]
						Y[1] = ∑Infilt_Obs[iSoil,infiltOutput.iT_TransSteady_Data[iSoil]]
						Plots.plot!(X, Y, seriestype=:scatter, label="TransSteady", color= :cyan, shape= :circle, markersize=10) 

						xlabel!(L"Time [s]")
						ylabel!(L"\sum infiltration \ [mm]")      
						
					Path = path.Plots_∑infilt_Tinfilt * "INFIL_" * option.infilt.Model * "_" * option.infilt.OutputDimension *  "_" * string(Id_Select[iSoil]) *  ".svg"

				Plots.savefig(Plot_∑infilt_Obs, Path)

			end # for iSoil=1:N_SoilSelect
			
			return
		end # PLOT_∑INFILT

end  # module plot
# ............................................................