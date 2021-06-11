# =============================================================
#		MODULE: plot
#
# =============================================================
module plot

	# =============================================================
	#		MODULE: lab
	# =============================================================
	module lab
		import ..cst, ..kunsat, ..param, ..wrc

		# using GLMakie
		using CairoMakie
	
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYDROPARAM
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, Id_Select, hydro, KunsatModel_Lab, Path_Plot_θΨK, Path_Model_Name; N_Se=1000)
			println("  ==  START: Plotting HydroParam  ==")
	
					
			# ===================== DATA =====================
			θ_Sim       = Array{Float64}(undef, (N_Se))
			θ_OtherData = Array{Float64}(undef, (N_Se))
			Kunsat_Sim  = Array{Float64}(undef, (N_Se))

			Ψ_θΨ_Min = 0.0

			for iZ = param.globalparam.N_iZ_Plot_Start: min(param.globalparam.N_iZ_Plot_End, N_SoilSelect)	
				Ψ_θΨ_Max = maximum(Ψ_θΨ[iZ,N_θΨ[iZ]]) + 100000.0

				Ψ_Sim = expm1.(range(log1p(Ψ_θΨ_Min), stop=log1p(Ψ_θΨ_Max), length=N_Se)) 

				θ_θΨ_Max = hydro.Φ[iZ]

				# Simulated 
					for iΨ = 1:N_Se
						θ_Sim[iΨ] = wrc.Ψ_2_θDual(optionₘ,Ψ_Sim[iΨ], iZ, hydro)

						Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(optionₘ, Ψ_Sim[iΨ], iZ, hydro)
					end # iΨ

				# == Title == 
					Title =  string(Id_Select[iZ]) * "_" * string(option.hydro.HydroModel) * "_" * string(option.hydro.σ_2_Ψm)

				# == Ticks ==
					Ticks = Int64.(Ψ_θΨ[iZ,1:N_θΨ[iZ]])
	
				#  ===================== PLOTTING =====================
					# CairoMakie.activate!(type = "svg")
					# inline!(true)
					Fig = Figure()
									
				#  == Plot_θ_Ψ  ==
					# Plot_θ_Ψ: General attributes
					# Axis1.xticks= (log1p.(cst.Mm_2_kPa * Ψ_θΨ[iZ,1:N_θΨ[iZ]])), Int64.(cst.Mm_2_kPa * Ψ_θΨ[iZ,1:N_θΨ[iZ]])
						Axis1 = Axis(Fig[1,1], resolution = (1000, 600))

						xlims!(Axis1, log1p.(cst.Mm_2_kPa * Ψ_θΨ_Min), log1p.(cst.Mm_2_kPa * Ψ_θΨ_Max * 1.1))

						# ylims!(Axis1, 0.0, θ_θΨ_Max)
						Axis1.xlabel = "ln(1 + Ψ) [kPa]"
						Axis1.ylabel =  "ln(1 + Ψ) [kPa]"
						Axis1.title = Title

						# Plot_θ_Ψ: Observed
						scatter!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_θΨ[iZ,1:N_θΨ[iZ]]), Float64.(θ_θΨ[iZ,1:N_θΨ[iZ]]), color=:red, markersize=10, marker = '■', label="Obs")

					# Plot_θ_Ψ: Simulated
							lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_Sim[1:N_Se], color=:blue, linewidth=2, label="Sim")

					# Plot_θ_Ψ: Total porosity point
						X = zeros(Float64,1)
						X[1] = 0.0
						Y = zeros(Float64,1)
						Y[1] = hydro.Φ[iZ]
						Label = 
						scatter!(Fig[1,1], log1p.(cst.Mm_2_kPa .* X), Y, color=:yellow, markersize=15, marker ="●", label="Φ")

					axislegend()

				# == Plot_K_Ψ  ==
				if option.hydro.KunsatΨ
					Axis2 = Axis(Fig[1,2])

					K_Ψ_Max = maximum(K_KΨ[iZ,1:N_KΨ[iZ]])
					xlims!(Axis2, log1p.(cst.Mm_2_kPa*Ψ_θΨ_Min), log1p.(Ψ_θΨ_Max*cst.Mm_2_kPa))
					ylims!(Axis2,  (log1p(0.0), log1p(hydro.Ks[iZ]* cst.MmS_2_CmH * 1.1)))
					Axis2.xlabel = "ln(1 + Ψ) [kPa]"
					Axis2.ylabel =  "ln ( 1 + K (Ψ) ) [cm h⁻¹]"              
	
					# Plot_K_Ψ: Ks                    
						X = zeros(Float64,1)
						X[1] = 0.0 
						Y = zeros(Float64,1)
						Y[1] = hydro.Ks[iZ]
						Label = "Ks_Max"
						scatter!(Fig[1,2], log1p.(X.* cst.Mm_2_kPa) , log1p.(Y.*cst.MmS_2_CmH) ,Y, color=:yellow, markersize=15, marker = '■', label=Label )

					# PlotK_Ψ: K(Ψ) obs
						X = Ψ_KΨ[iZ,1:N_KΨ[iZ]]
						Y = K_KΨ[iZ,1:N_KΨ[iZ]]
						Label = "Obs"
						scatter!(Fig[1,2], log1p.(X.* cst.Mm_2_kPa) , log1p.(Y.*cst.MmS_2_CmH) ,Y, color=:red, markersize=10, marker = '■', label=Label )

					# Plot_K_Ψ: K(Ψ) sim
						X = Ψ_Sim
						Y = Kunsat_Sim 
						Label = "Sim"
						lines!(Fig[1,2], log1p.(Ψ_Sim.*cst.Mm_2_kPa), log1p.(Kunsat_Sim.*cst.MmS_2_CmH), color=:blue, label=Label)

						axislegend()

					# General attributes
					
					#   Axis2.xticks(log1p.(Ψ_Sim.*cst.Mm_2_kPa))
					#   , Int64.(cst.Mm_2_kPa * Ψ_θΨ[iZ,1:N_θΨ[iZ]])
					
					# Plot1 = Plots.plot(Plot1,Plot_Θψ, Plot_kΘ)
				# else
				#    Plot1 = Plots.plot(Plot1, Plot_Θψ)

					
				end # option.hydro.KunsatΨ 

				# Fig[3, 1] = Legend(Fig, Axis1, "PLOTS", orientation=:horizontal)

				# Path = path.plotSoilwater.Plot_θΨK * "Lab_ThetaH_" * Title * ".svg" 
				
				Path = Path_Plot_θΨK * "Lab_ThetaH_" * string(Path_Model_Name) * "_" * string(Id_Select[iZ]) * ".svg" 
	
				save(Path, Fig)
	
				display(Fig)

				println("    ~  $(Path) ~")
			end # for iZ

			println("  ==  END: Plotting HydroParam  == \n")		
		return nothing

	end  # module lab
	# ............................................................



	# =============================================================
	#		MODULE: psd
	# =============================================================
	module psd
		using Plots, Plots.PlotMeasures, LaTeXStrings
		import ...wrc, ...kunsat, ...cst, ...param, ...psdThetar, ...psdFunc, ...bestFunc
		export PLOT_θr, PLOT_IMP_MODEL, PLOT_PSD_θΨ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_θr
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_θr(∑Psd, N_SoilSelect, hydro, hydroPsd, Path)
			println("  ==  START: Plotting PLOT_θr  ==")

			# Sorting ascending order with clay fraction
				Array      = zeros(Float64, 3, length(∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size]))
				Array      = zeros(Float64, (3, N_SoilSelect))
			
				Array[1,:] = ∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size] # Clay fraction
				Array[2,:] = hydroPsd.θr[1:N_SoilSelect]
				Array[3,:] = hydro.θr[1:N_SoilSelect]
				Array      = sortslices(Array, dims=2)
				Clay       = Array[1,:] # Clay fraction
				θr_Psd     = Array[2,:]
				θr         = Array[3,:]
			
			# Minimum and maximum value
				θr_Min = 0.01 

				θr_Max = maximum(hydroPsd.θr_Max) + 0.05
				Clay_Min = 0.1
				Clay_Max = maximum(∑Psd[1:N_SoilSelect, param.psd.Psd_2_θr_Size]) + 0.05
			
			# PLOT 1 <>=<>=<>=<>=<>=<>
			pgfplotsx()
				# Plot θr(Clay)
					X = Clay
					Y = θr
					Plot_θr = Plots.plot(X, Y, seriestype=:scatter, label=L"\theta _{r}", color= :violet, shape= :square, markersize=4, legend=:bottomright, size=(5000,400))

				# Plot θr_psd(Clay)
					X = Clay
					Y = θr_Psd
					Plots.plot!(X ,Y, seriestype=:line, label=L"\theta _{r psd}", color= :blue, lw=2)
			
				# General attributes
					xlabel!(L"Clay \ [g \ g^{-1}]")                         
					ylabel!(L"\theta _{r} [cm^{3} \ cm^{-3}]")
					Plots.plot!(xlims= (Clay_Min, Clay_Max), ylims= (θr_Min, θr_Max))

			# PLOT 2 <>=<>=<>=<>=<>=<>
				# Plot θr_Psd(θr)
					X = θr
					Y = θr_Psd
					Plot_θr_Psd = Plots.plot(X ,Y, seriestype=:scatter, color=:violet, shape=:square, markersize=4, size=(800,400))
					
				# 1:1 line
					X = range(θr_Min, stop=θr_Max, length=10)
					Y = X
					Label = "1:1"
					Plots.plot!(X, Y, seriestype=:line, label= Label , color= :black, linestyle= :dot, lw=2)

				# General attributes
					xlabel!(L"\theta _{r} [cm^3 cm^{-3}]")
					ylabel!(L"\theta _{r \ psd} [cm^3 cm^{-3}]")
					Plots.plot!(xlims= (θr_Min, θr_Max), ylims= (θr_Min, θr_Max))

			Plot = Plots.plot(Plot_θr, Plot_θr_Psd)
			Plots.savefig(Plot, Path)
			println("    ~  $(Path) ~")
		
		println("  ==  END: Plotting PLOT_θr  == \n")
		return
		end # function: PLOT_θr


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_IMP_MODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_IMP_MODEL(Id_Select, Rpart, N_Psd, ∑Psd, Psd, N_SoilSelect, hydro, paramPsd, Path)
			println("  ==  START: PLOT_IMP_MODEL  ==")	

			for iZ = param.globalparam.N_iZ_Plot_Start: param.globalparam.N_iZ_Plot_End
				Rpart_Min = minimum(Rpart[iZ,1:N_Psd[iZ]])
				Rpart_Max = maximum(Rpart[iZ,1:N_Psd[iZ]])

				∑Psd_Min  = minimum(∑Psd[iZ,1:N_Psd[iZ]])
				∑Psd_Max  = maximum(∑Psd[iZ,1:N_Psd[iZ]])

				Psd_Min  = minimum(Psd[iZ,1:N_Psd[iZ]])
				Psd_Max  = maximum(Psd[iZ,1:N_Psd[iZ]])

				IntergranularMixing = zeros(Float64, N_Psd[iZ])
				ξ = zeros(Float64, N_Psd[iZ])
				for iRpart = 1:N_Psd[iZ]
					ξ[iRpart] = psdFunc.imp.INTERGRANULARMIXING(Rpart[iZ,iRpart], paramPsd.ξ1[iZ], paramPsd.ξ2[iZ])

					IntergranularMixing[iRpart] = (Rpart[iZ, iRpart] ^ -ξ[iRpart]) 
				end # for iRpart = 1:N_Psd[iZ]

				# << PLOT 1 >>
					# Plot_∑Psd_Rpart
						X = Rpart[iZ,1:N_Psd[iZ]]
						Y = ∑Psd[iZ,1:N_Psd[iZ]]
						Plot_∑Psd_Rpart = Plots.plot(X ,Y, seriestype=:scatter, color= :teal, shape= :square, markersize= 4, size=(800,400))
						Plots.plot!(X ,Y, seriestype=:line, color= :teal)

					# Plot_∑Psd_Rpart: General attributes
						xlabel!(L"R_{part} [mm]")
						ylabel!(L"\sum \ PSD")
						Plots.plot!(xlims= (Rpart_Min, Rpart_Max), ylims= (∑Psd_Min, ∑Psd_Max), xscale= :log10)

				# << PLOT 2 >>
					# Plot_Psd_Rpart
						X = Rpart[iZ,1:N_Psd[iZ]]
						Y = Psd[iZ,1:N_Psd[iZ]]
						# Plot_Rpart_Psd = Plots.plot(X ,Y, seriestype=:scatter, color= :blue, shape= :square, markersize= 4, size=(800,400))
						Plot_Rpart_Psd = Plots.plot(X ,Y, seriestype=:line, color= :blue, shape= :square, markersize= 4, size=(800,400))

						xlabel!(L"R_{part} \ [mm]")
						ylabel!(L"PSD [mm]")
						Plots.plot!(xlims= (Rpart_Min, Rpart_Max), ylims= (Psd_Min, Psd_Max), xscale= :log10)


				# << PLOT 3 >>
					# Plot NormMixing_Rpart
						X = Rpart[iZ,1:N_Psd[iZ]]
						Y = IntergranularMixing[1:N_Psd[iZ]] / maximum( IntergranularMixing[1:N_Psd[iZ]] )
						Plot_NormMixing_Rpart = Plots.plot(X, Y, seriestype=:line, color= :green)

						xlabel!(L"R_{part} \ [mm]")
						ylabel!(L"R_{part}^{-\xi(R_{Part})}")
						Plots.plot!(xlims= (Rpart_Min, Rpart_Max), ylims= (0.0, 1.1), xscale= :log10)


				Plot = Plots.plot(Plot_∑Psd_Rpart, Plot_Rpart_Psd, Plot_NormMixing_Rpart, layout = (3,1))
				Path = Path * "IMP_" * string(option.hydro.HydroModel) * "_" *string(Id_Select[iZ]) * ".svg"
				Plots.savefig(Plot, Path)
				println("    ~  $(Path) ~")
			end # for iZ
			println("  ==  END: PLOT_IMP_MODEL  == \n")
			return	
		end # function: PLOT_IMP_MODEL


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_IMP_ΘΨ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_PSD_θΨ(Ψ_θΨ, Ψ_Rpart, θ_θΨ, θ_Rpart, N_θΨ, N_SoilSelect, N_Psd, Id_Select, hydroPsd, hydro, Path; N_Se= 100)
		
			println("  ==  START: Plotting PLOT_PSD_θΨ  ==")

			θ_θΨ_Psd = Array{Float64}(undef, (N_Se))

			for iZ = param.globalparam.N_iZ_Plot_Start:min(param.globalparam.N_iZ_Plot_End, N_SoilSelect)
				# Range of plots
					Ψ_θΨ_Min = 10.0 ^ 0.0 # [mm]

					Ψ_θΨ_Max = 150000 + 100000 # [mm]

					Ψ_Sim = 10.0 .^ range(log(Ψ_θΨ_Min), stop=log(Ψ_θΨ_Max), length=N_Se)

					θ_θΨ_Max = hydroPsd.Φ[iZ] + 0.1

				# Simulated 
					for iΨ = 1:N_Se
						θ_θΨ_Psd[iΨ] = wrc.Ψ_2_θDual(optionₘ,Ψ_Sim[iΨ], iZ, hydroPsd)
					end # iΨ 

				# Plot_θ_Ψ: Psd model fitting for e.g. Kosugi model
					X = Ψ_Sim[1:N_Se] .* cst.Mm_2_Cm
					Y = θ_θΨ_Psd[1:N_Se]
					Label = "PsdKosugi"
					Plot_θ_Ψ_Psd = Plots.plot(X ,Y, seriestype=:line, label=Label, color= :blue, lw=2)

				# Plot_θ_Ψ: Psd model points
					X = Ψ_Rpart[iZ, 1:N_Psd[iZ]] .* cst.Mm_2_Cm
					Y = θ_Rpart[iZ, 1:N_Psd[iZ]]
					Label = "PsdModel"
					Plot_θ_Ψ_Psd = Plots.plot!(X ,Y, seriestype=:scatter, label=Label, color= :violet, shape= :circle, markersize=4)

				# Plot_θ_Ψ: Observed
				if option.run.HydroLabθΨ ≠ :No 
					X = Ψ_θΨ[iZ,1:N_θΨ[iZ]] .* cst.Mm_2_Cm
					Y = θ_θΨ[iZ,1:N_θΨ[iZ]]
					Label = "LabObs"
					Plot_θ_Ψ = Plots.plot!(X ,Y, seriestype=:scatter, label=Label, color= :red, shape= :square, markersize=4)

					# Plot_θ_Ψ: Total porosity point
						X = zeros(Float64,1)
						X[1] = Ψ_θΨ_Min * cst.Mm_2_Cm
						Y = zeros(Float64,1)
						Y[1] = hydro.Φ[iZ]
						Label = "\$ \\phi \$"
						Plots.plot!(X, Y, seriestype=:scatter, label= Label, color= :green, shape= :square, markersize=4) 
				end


				# Plot_θ_Ψ: General attributes
					xlabel!(L"\psi \ [cm]")
					ylabel!(L"\theta \ [cm^3 cm^{-3}]")
					Plots.plot!(xlims =(Ψ_θΨ_Min*cst.Mm_2_Cm, Ψ_θΨ_Max*cst.Mm_2_Cm), ylims =(0.0, θ_θΨ_Max), xscale= :log10, size=(800,400))

				Path = Path * "Psd_ThetaH_" * string(option.hydro.HydroModel) * "_" *string(Id_Select[iZ]) * ".svg"     
				Plot = Plots.plot(Plot_θ_Ψ_Psd)
				Plots.savefig(Plot, Path)
				println("    ~  $(Path) ~")
			end # iZ
		println("  ==  END: Plotting PLOT_PSD_θΨ  == \n")
		return	nothing	
		end # function PLOT_IMP_ΘΨ

	end  # module: psd
	# ............................................................


	# =============================================================
	#		MODULE: infilt
	# =============================================================
	module infilt
		import ...wrc, ...kunsat, ...cst, ...param, ...psdThetar, ...psdFunc, ...bestFunc, ...sorptivity
		using Plots, Plots.PlotMeasures, LaTeXStrings
		export  PLOT_∑INFILT, PPLOT_∑INFILT_θΨ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_∑INFILT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_∑INFILT(Id_Select, N_Infilt, N_SoilSelect, ∑Infilt_Obs, Tinfilt, ∑Infilt_3D, ∑Infilt_1D, infiltOutput, Path)
		println("  ==  START: PLOT_∑INFILT  == \n")
		
			for iZ = param.globalparam.N_iZ_Plot_Start: min(param.globalparam.N_iZ_Plot_End, N_SoilSelect)	
				# << PLOT 1 >>
					Title = " iZ= $(Id_Select[iZ])"
					# Plot_∑infilt_Obs

						Label ="Obs_$(string(option.infilt.DataSingleDoubleRing))_Ring"
						X = Tinfilt[iZ,1:N_Infilt[iZ]] / 60.0
						Y = ∑Infilt_Obs[iZ,1:N_Infilt[iZ]]
						Plot_∑infilt_Obs = Plots.plot(X, Y, seriestype=:scatter, label=Label, color= :red, shape= :square, markersize=4, marker = (Plots.stroke(1, :red))) 

					# Plot_∑infilt_Sim
						Label = "Sim_3D"
						X = Tinfilt[iZ,1:N_Infilt[iZ]] / 60.0
						Y = ∑Infilt_3D[iZ,1:N_Infilt[iZ]]
						Plots.plot!(X, Y, seriestype=:line, label=Label, color= :blue, shape= :square, markersize=4, marker = (Plots.stroke(1, :blue))) 


						Label = "Sim_1D"
						Y2 = ∑Infilt_1D[iZ,1:N_Infilt[iZ]]
						Plots.plot!(X, Y2, seriestype=:line, label=Label, color= :green, shape= :square, markersize=4,  marker = (Plots.stroke(1, :green))) 

					# TransSteady
						Label="T_TransSteady"
						X = zeros(Float64,1)
						Y = zeros(Float64,1)
						X[1] = Tinfilt[iZ,infiltOutput.iT_TransSteady_Data[iZ]] / 60.0
						Y[1] = ∑Infilt_Obs[iZ,infiltOutput.iT_TransSteady_Data[iZ]]
						Plots.plot!(X, Y, seriestype=:scatter, label=Label, color= :violet, shape= :circle, markersize=10, title=Title) 

						Plots.xlabel!(L"Time [minutes]")
						Plots.ylabel!(L"\sum infiltration \ [mm]")      
						
					Path = Path * "INFIL_" * string(option.infilt.Model)  *  "_" * string(Id_Select[iZ]) *  ".svg"

				Plots.savefig(Plot_∑infilt_Obs, Path)
				println("    ~  $(Path) ~")

			end # for iZ=1:N_SoilSelect
		println("  ==  END: PLOT_∑INFILT  == \n")
		return nothing
		end # PLOT_∑INFILT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_∑INFILT_θΨ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_∑INFILT_θΨ(hydroInfilt, Id_Select, N_SoilSelect, Param, Path; hydro=[], N_Se=100)
		println("  ==  START: PLOT_∑INFILT_θΨ  ==")

			θ_Infilt      = Array{Float64}(undef, (N_Se))
			θ_Obs         = Array{Float64}(undef, (N_Se))
			Kunsat_Infilt = Array{Float64}(undef, (N_Se))
			Kunsat_Obs    = Array{Float64}(undef, (N_Se))

			for iZ = param.globalparam.N_iZ_Plot_Start: param.globalparam.N_iZ_Plot_End	
				Ψ_θΨ_Min = 10.0 ^ -2 # [mm]

				Ψ_θΨ_Max = 200000.0 * 10.0 # [mm]

				Ψ = 10.0 .^ range(log(Ψ_θΨ_Min), stop=log(Ψ_θΨ_Max), length=N_Se)

				θ_θΨ_Max = hydroInfilt.Φ[iZ] + 0.1

				if option.run.HydroLabθΨ ≠ :No && option.hydro.KunsatΨ
					K_Ψ_Max = max(hydroInfilt.Ks[iZ], hydro.Ks[iZ]) * 1.1
				else
					K_Ψ_Max = hydroInfilt.Ks[iZ] * 1.1
				end #  option.hydro.KunsatΨ

				for iΨ = 1:N_Se
					θ_Infilt[iΨ] = wrc.Ψ_2_θDual(optionₘ,Ψ[iΨ], iZ, hydroInfilt)

					Kunsat_Infilt[iΨ] = kunsat.Ψ_2_KUNSAT(optionₘ, Ψ[iΨ], iZ, hydroInfilt)

					if option.run.HydroLabθΨ ≠ :No
						θ_Obs[iΨ] = wrc.Ψ_2_θDual(optionₘ,Ψ[iΨ], iZ, hydro)

						if option.run.HydroLabθΨ ≠ :No && option.hydro.KunsatΨ
							Kunsat_Obs[iΨ] = kunsat.Ψ_2_KUNSAT(optionₘ, Ψ[iΨ], iZ, hydro)
						end # option.hydro.KunsatΨ		
					end # option.run.HydroLabθΨ ≠ :No

				end # iΨ 

				#PLOT 1:  Plot_θ_Ψ
					# Plot_θ_Ψ: Simulated Infiltration
						X = Ψ[1:N_Se] .* cst.Mm_2_Cm
						Y = θ_Infilt[1:N_Se]
						Label = "Infilt"
						Plots.plot(X, Y, seriestype=:line, label=Label, color= :blue, lw=2)

					# Plot_θ_Ψ: Observed
					if option.run.HydroLabθΨ ≠ :No
						X = Ψ[1:N_Se] .* cst.Mm_2_Cm
						Y = θ_Obs[1:N_Se]
						Label = "Obs"
						Plot_θ_Ψ = Plots.plot!(X ,Y, seriestype=:line, label=Label, color= :red, lw=2)
					end # option.run.HydroLabθΨ ≠ :No

					# Plot_θ_Ψ: General attributes
						Plots.xlabel!("\\psi [cm]")
						Plots.ylabel!(L"\theta \ [cm^3 cm^{-3}]")
						Plots.plot!(xlims =(10.0*Ψ_θΨ_Min*cst.Mm_2_Cm, Ψ_θΨ_Max*cst.Mm_2_Cm), ylims =(0.0, θ_θΨ_Max), xscale= :log10, size=(800,400), legend=:bottomleft)

					# PLOT2: Kunsat
						# Plot_K_Ψ: Obs K_Ψ
						X = Ψ[1:N_Se] .* cst.Mm_2_Cm
						Y = Kunsat_Infilt[1:N_Se] .* cst.MmS_2_CmH
						Label = "Infilt"
						Plot_K_Ψ = Plots.plot(X, Y, seriestype=:line, label=Label, color= :blue, lw=2)

						# Plot_K_Ψ: Sim K_Ψ
						if option.run.HydroLabθΨ ≠ :No && option.hydro.KunsatΨ
							X = Ψ[1:N_Se] .* cst.Mm_2_Cm
							Y = Kunsat_Obs[1:N_Se] .* cst.MmS_2_CmH
							Label = "Obs"
							Plots.plot!(X, Y, seriestype=:line, label=Label, color= :red, lw=2)
						end # option.hydro.KunsatΨ

						# General attributes
							Plots.xlabel!("\\psi [cm]")
							Plots.ylabel!(L" K (\psi) \ [cm \ h^{-1}]")
							Plots.plot!(xlims = (Ψ_θΨ_Min*cst.Mm_2_Cm, Ψ_θΨ_Max*cst.Mm_2_Cm), ylims = (10^-2.0, K_Ψ_Max * cst.MmS_2_CmH), xscale= :log10,  yscale= :log10, legend=:bottomleft, size=(800,400))

					Path = Path * "Infilt_ThetaH_" * string(option.hydro.HydroModel) * "_" *string(Id_Select[iZ]) * ".svg"     
					Plot = Plots.plot(Plot_θ_Ψ, Plot_K_Ψ)
					Plots.savefig(Plot, Path)
					println("    ~  $(Path) ~")
			end # iZ

		println("  ==  END: PLOT_∑INFILT_θΨ  == \n")
		return
		end  # function: PLOT_∑INFILT_θΨ


end  # module plot
# ............................................................
