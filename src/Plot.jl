# =============================================================
#		MODULE: plot
#
# =============================================================
module plot

	# =============================================================
	#		MODULE: lab
	# =============================================================
	module lab
		import ...cst, ...kunsat, ...wrc
		using CairoMakie

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYDROPARAM
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYDROPARAM(hydro, hydroOther, IdSelect, K_KΨobs, NiZ, N_KΨobs, N_θΨobs, optim, option, param, path, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs; N_Se=1000)

				println("  ==  START: Plotting HydroParam  ==")

				# ===================== DATA =====================
				θ_Sim             = fill(0.0,N_Se)
				Kunsat_Sim        = fill(0.0,N_Se)

				Ψ_θΨobs_Min = 0.0
				for iZ = param.globalparam.N_iZ_Plot_Start:param.globalparam.N_iZ_Plot_End
					Ψ_θΨobs_Max = maximum(Ψ_θΨobs[iZ,N_θΨobs[iZ]]) + 100000.0

					Ψ_Sim = expm1.(range(log1p(Ψ_θΨobs_Min), stop=log1p(Ψ_θΨobs_Max), length=N_Se)) 

					θ_θΨobs_Max = hydro.Φ[iZ]

					# Simulated 
						for iΨ = 1:N_Se
							θ_Sim[iΨ] = wrc.Ψ_2_θDual(option.hydro, Ψ_Sim[iΨ], iZ, hydro)
							Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(option.hydro, Ψ_Sim[iΨ], iZ, hydro)
						end # iΨ = 1:N_Se

					# _______________________ START: Plotting _______________________
								
					Fig = Figure(backgroundcolor=RGBf0(0.98, 0.98, 0.98), resolution = (2500, 1000),  font="Sans", fontsize=16)

					Title = "iZ= $(IdSelect[iZ]) " * "θ(Ψ) Nse_θΨ=" * string(round(hydroOther.Nse_θΨ[iZ], digits=2)) * "; Nse_KΨ=" * string(round(hydroOther.Nse_KΨ[iZ], digits=2)) * "; Wilmot_θΨ=" *  string(round(hydroOther.NseWilmot_θΨ[iZ],digits=2)) * "; Wilmot_KΨ=" * string(round(hydroOther.NseWilmot_KΨ[iZ], digits=2))
					
					#  == Plot_θ_Ψ  ==
						Axis1 = Axis(Fig[1,1], title=Title, titlesize=24, xlabel="ln(1 + Ψ) [kPa]", ylabel="θ [mm³ mm⁻³]", xlabelsize=10, backgroundcolor=:white)

						xlims!(Axis1, log1p.(cst.Mm_2_kPa * Ψ_θΨobs_Min), log1p.(cst.Mm_2_kPa * Ψ_θΨobs_Max ))
						ylims!(Axis1, 0.0, 0.75)
						Axis1.xticks = (log1p.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]]), string.( floor.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]], digits=1)))

						Fig_θΨobs = scatter!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_θΨobs[iZ,1:N_θΨobs[iZ]]), Float64.(θ_θΨobs[iZ,1:N_θΨobs[iZ]]), color=:red, markersize=25, marker = '■')

						Fig_θΨsim = lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_Sim[1:N_Se], color=:blue, linewidth=3)

					# Plot_θ_Ψ: Total porosity point
						Fig_TotalPorosity = scatter!(Fig[1,1], [log1p.(cst.Mm_2_kPa .* 0.0)], [hydro.Φ[iZ]], color=:green, markersize=25, marker ='●')

					# == Plot_K_Ψ  ==
					# If Ks is not computed it is computed from Ks_Model

						Axis2 = Axis(Fig[1,2], title="K(Ψ)", titlesize=24, xlabel = "ln(1 + Ψ) [kPa]", ylabel = "ln (1 + K (Ψ)) [mm h⁻¹]")

						xlims!(Axis2, log1p.(cst.Mm_2_kPa * Ψ_θΨobs_Min), log1p.(cst.Mm_2_kPa * Ψ_θΨobs_Max))

						ylims!(Axis2, 0.0, log1p(hydro.Ks[iZ]*cst.MmS_2_MmH))
						
						Axis2.xticks = (log1p.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]]), string.(floor.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]], digits=1)))
						Yticks = 1:1:6
						Axis2.yticks = (Yticks,string.(Yticks))

						if option.data.Kθ
							Fig_Kθobs = scatter!(Fig[1,2], log1p.(Ψ_KΨobs[iZ,1:N_KΨobs[iZ]].*cst.Mm_2_kPa), log1p.(K_KΨobs[iZ,1:N_KΨobs[iZ]].*cst.MmS_2_MmH), color=:red, markersize=25, marker = '■')
						end

						Fig_Kθsim = lines!(Fig[1,2], log1p.(Ψ_Sim[1:N_Se].*cst.Mm_2_kPa), log1p.(Kunsat_Sim[1:N_Se] .* cst.MmS_2_MmH), color=:blue, linewidth=3)

						Fig_Ks = scatter!(Fig[1,2], [log1p.(cst.Mm_2_kPa .* 0.0)], [log1p(hydro.Ks[iZ] * cst.cst.MmS_2_MmH)], color=:green, markersize=25, marker ='●')


						if option.data.Kθ
							Leg = Fig[1, end+1] = Legend(Fig, [Fig_θΨobs, Fig_θΨsim, Fig_TotalPorosity, Fig_Kθobs, Fig_Kθsim, Fig_Ks], ["θobs(Ψ)", "θsim(Ψ)", "Φ", "Kobs(Ψ)", "Ksim(Ψ)", "Ksₛᵢₘ"])
						else
							Leg = Fig[1, end+1] = Legend(Fig, [Fig_θΨobs, Fig_θΨsim, Fig_TotalPorosity, Fig_Kθsim, Fig_Ks], ["θobs(Ψ)", "θsim(Ψ)", "Φ", "Ksim(Ψ)", "Ksₛᵢₘ"])
						end

					Fig[2, 1:2] = Leg
					trim!(Fig.layout)
					Leg.orientation = :horizontal
					Leg.tellheight = true
					
					Path = path.plotSoilwater.Plot_θΨK * "Lab_ThetaH_" * string(path.option.ModelName) * "_" * string(IdSelect[iZ]) * ".svg" 
					save(Path, Fig)
	
					# Displaying figure in VScode
					# if option.other.PlotVscode
					# 	display(Fig)
					# end

					# println("    ~  $(Path) ~")
				
				end # for iZ
				
			# ------------------------END: Plotting---------------------------  
			println("  ==  END: Plotting HydroParam  == \n")		
			return nothing
			end  # function: HYDROPARAM
	
	end  # module lab
	# ............................................................


	# =============================================================
	#		module: ksmodel

	# =============================================================
	module ksmodel
		import ...cst
		using CairoMakie
		# using GLMakie
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KSMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSMODEL(hydro, KₛModel, Ksₒᵦₛ, NameSim::String, Path::String, θrₒᵦₛ, θsMacMatₒᵦₛ, σₒᵦₛ)
			Ks_Min = minimum([minimum(Ksₒᵦₛ), minimum(KₛModel)])
			Ks_Max = maximum([maximum(Ksₒᵦₛ), maximum(KₛModel)])

			Ks_Max = 0.099371778 # mm/s
			
			CairoMakie.activate!()

			Fig = Figure(backgroundcolor=RGBf0(0.98, 0.98, 0.98), resolution = (2500, 1000),  font="Sans", fontsize=20, xgridstyle=:dash, ygridstyle=:dash, xtickalign=1, ytickalign=1)
					
			#  == Plot_θ_Ψ  == 
				Axis1 = Axis(Fig[1,1], title="KsModel_" * NameSim, titlesize=25, xlabel="ln (1 + Ks_Obs) [mm hour⁻¹ ]", ylabel="ln (1 + KsModel_Sim) [mm hour⁻¹ ]", xlabelsize=25,  ylabelsize=25, xticksize=20,  yticksize=20)

				xlims!(Axis1, log1p.(0.0), log1p.(Ks_Max * cst.MmS_2_MmH))
				ylims!(Axis1, log1p.(0.0), log1p.(Ks_Max * cst.MmS_2_MmH))

				KsTicks = expm1.(range(log1p(0.0), stop=log1p(Ks_Max * cst.MmS_2_MmH), length=10)) 
				Axis1.xticks = (log1p.(KsTicks), string.( floor.(KsTicks, digits=1)))
				Axis1.yticks = (log1p.(KsTicks), string.( floor.(KsTicks, digits=1)))

				ΔΘsMacΘr = θsMacMatₒᵦₛ .-  θrₒᵦₛ

				Fig_Ks = scatter!(Fig[1,1], log1p.(Ksₒᵦₛ .* cst.MmS_2_MmH), log1p.(KₛModel .* cst.MmS_2_MmH), color=σₒᵦₛ, markersize=125.0*ΔΘsMacΘr, marker =:circle)

				Colorbar(Fig[1, 2], limits=(minimum(σₒᵦₛ), maximum(σₒᵦₛ)), colormap =:viridis, label="σ[-]", vertical = true)

				Line = range(log1p(Ks_Min.* cst.MmS_2_MmH), stop=log1p(Ks_Max.* cst.MmS_2_MmH), length=10) 

				Fig_Ks = lines!(Fig[1,1], Line, Line, color=:blue)

				# Leg1 = Colorbar(Fig, Fig_Ks, label = "Theta", ticklabelsize = 14, labelpadding = 5, width = 10)
				
				trim!(Fig.layout)

			Fig[1, 1] = Axis1
   		# Fig[1, 2] = Leg1

			Pathₛ = Path * "_" * NameSim * ".svg" 

			save(Pathₛ, Fig)
			display(Fig)

		return nothing
		end  # function: KSMODEL
		# ------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KsModel_3D
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KsModel_3D(hydro, NiZ, path)

			GLMakie.activate!()

			Fig = Figure(backgroundcolor=RGBf0(0.98, 0.98, 0.98), resolution = (2500, 1000),  font="Sans", fontsize=35)

			Ks_3D = surface(Fig[1,1], hydro.θs[1:NiZ] .- hydro.θr[1:NiZ], hydro.σ[1:NiZ], hydro.Ks[1:NiZ]; shading=false, colormap = :deep, axis = (show_axis = false,))

			Path = path.plotSoilwater.Plot_θΨK * "//KsModel//" 
			mkpath(Path)
			Path = Path * "KsModel_3D.svg" 
			save(Path, Fig)
			display(Fig)

		return nothing
		end  # function: KsModel_3D
		# ------------------------------------------------------------------
		
	end  # module: ksmodel

	# ............................................................

	# =============================================================
	#		MODULE: psd
	# =============================================================
	module psd
		using Plots, Plots.PlotMeasures, LaTeXStrings
		import ...wrc, ...kunsat, ...cst, ...psdThetar, ...psdFunc, ...bestFunc
		export PLOT_θr, PLOT_IMP_MODEL, PLOT_PSD_θΨ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_θr
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PLOT_θr(∑Psd, hydro, hydroPsd, NiZ, param, Path)
				println("  ==  START: Plotting PLOT_θr  ==")

				# Sorting ascending order with clay fraction
					Array      = zeros(Float64, 3, length(∑Psd[1:NiZ, param.psd.Psd_2_θr_Size]))
					Array      = zeros(Float64, (3, NiZ))
				
					Array[1,:] = ∑Psd[1:NiZ, param.psd.Psd_2_θr_Size] # Clay fraction
					Array[2,:] = hydroPsd.θr[1:NiZ]
					Array[3,:] = hydro.θr[1:NiZ]
					Array      = sortslices(Array, dims=2)
					Clay       = Array[1,:] # Clay fraction
					θr_Psd     = Array[2,:]
					θr         = Array[3,:]
				
				# Minimum and maximum value
					θr_Min = 0.01 

					θr_Max = maximum(hydroPsd.θr_Max) + 0.05
					Clay_Min = 0.1
					Clay_Max = maximum(∑Psd[1:NiZ, param.psd.Psd_2_θr_Size]) + 0.05
				
				# PLOT 1 <>=<>=<>=<>=<>=<>
					# pgfplotsx()
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
			return nothing
			end # function: PLOT_θr


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_IMP_MODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PLOT_IMP_MODEL(∑Psd, hydro, IdSelect, NiZ, N_Psd, option, param, Path, Psd, Rpart)
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
						# added at a later stage
						ξ2 = psdFunc.imp.∑PSD_2_ξ2(∑Psd[param.psd.imp.∑Psd_2_ξ2_Size], param)

						ξ[iRpart] = psdFunc.imp.INTERGRANULARMIXING(param, Rpart[iZ,iRpart], param.psd.imp.ξ1, ξ2)

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
					Path₀ = Path * "IMP_" * string(option.hydro.HydroModel⍰) * "_" *string(IdSelect[iZ]) * ".svg"
					Plots.savefig(Plot, Path₀)
					println("    ~  $(Path₀) ~")
				end # for iZ
			println("  ==  END: PLOT_IMP_MODEL  == \n")
			return nothing	
			end # function: PLOT_IMP_MODEL


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_IMP_ΘΨ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PLOT_PSD_θΨ(hydro, hydroPsd, IdSelect, NiZ, N_Psd, N_θΨobs, option, param, Path, θ_Rpart, θ_θΨobs, Ψ_Rpart, Ψ_θΨobs; N_Se= 100)
			
				println("  ==  START: Plotting PLOT_PSD_θΨ  ==")

				θ_θΨobs_Psd = fill(0.0::Float64, (N_Se))

				for iZ = param.globalparam.N_iZ_Plot_Start:param.globalparam.N_iZ_Plot_End
					# Range of plots
						Ψ_θΨobs_Min = 10.0 ^ 0.0 # [mm]

						Ψ_θΨobs_Max = 150000 + 100000 # [mm]

						Ψ_Sim = 10.0 .^ range(log(Ψ_θΨobs_Min), stop=log(Ψ_θΨobs_Max), length=N_Se)

						θ_θΨobs_Max = hydroPsd.Φ[iZ] + 0.1

					# Simulated 
						for iΨ = 1:N_Se
							θ_θΨobs_Psd[iΨ] = wrc.Ψ_2_θDual(option.psd,Ψ_Sim[iΨ], iZ, hydroPsd)
						end # iΨ 

					# Plot_θ_Ψ: Psd model fitting for e.g. Kosugi model
						X = Ψ_Sim[1:N_Se] .* cst.Mm_2_Cm
						Y = θ_θΨobs_Psd[1:N_Se]
						Label = "PsdKosugi"
						Plot_θ_Ψ_Psd = Plots.plot(X ,Y, seriestype=:line, label=Label, color= :blue, lw=2)

					# Plot_θ_Ψ: Psd model points
						X = Ψ_Rpart[iZ, 1:N_Psd[iZ]] .* cst.Mm_2_Cm
						Y = θ_Rpart[iZ, 1:N_Psd[iZ]]
						Label = "PsdModel"
						Plot_θ_Ψ_Psd = Plots.plot!(X ,Y, seriestype=:scatter, label=Label, color= :violet, shape= :circle, markersize=4)

					# Plot_θ_Ψ: Observed
					if option.run.HydroLabθΨ⍰ ≠ "No" 
						X = Ψ_θΨobs[iZ,1:N_θΨobs[iZ]] .* cst.Mm_2_Cm
						Y = θ_θΨobs[iZ,1:N_θΨobs[iZ]]
						Label = "LabObs"
						Plot_θ_Ψ = Plots.plot!(X ,Y, seriestype=:scatter, label=Label, color= :red, shape= :square, markersize=4)

						# Plot_θ_Ψ: Total porosity point
							X = zeros(Float64,1)
							X[1] = Ψ_θΨobs_Min * cst.Mm_2_Cm
							Y = zeros(Float64,1)
							Y[1] = hydro.Φ[iZ]
							Label = "\$ \\phi \$"
							Plots.plot!(X, Y, seriestype=:scatter, label= Label, color= :green, shape= :square, markersize=4) 
					end

					# Plot_θ_Ψ: General attributes
						xlabel!(L"\psi \ [cm]")
						ylabel!(L"\theta \ [cm^3 cm^{-3}]")
						Plots.plot!(xlims =(Ψ_θΨobs_Min*cst.Mm_2_Cm, Ψ_θΨobs_Max*cst.Mm_2_Cm), ylims =(0.0, θ_θΨobs_Max), xscale= :log10, size=(800,400))

					Path₀ = Path * "Psd_ThetaH_" * string(option.hydro.HydroModel⍰) * "_" *string(IdSelect[iZ]) * ".svg"     
					Plot = Plots.plot(Plot_θ_Ψ_Psd)
					Plots.savefig(Plot, Path₀)
					println("    ~  $(Path₀) ~")
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
		import ...wrc, ...kunsat, ...cst, ...psdThetar, ...psdFunc, ...bestFunc, ...sorptivity
		using Plots, Plots.PlotMeasures, LaTeXStrings
		export  PLOT_∑INFILT, PLOT_∑INFILT_θΨ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_∑INFILT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PLOT_∑INFILT(∑Infilt_1D, ∑Infilt_3D, ∑Infilt_Obs, IdSelect, infiltOutput, N_Infilt, NiZ, option, param, Path, Tinfilt)
			println("  ==  START: PLOT_∑INFILT  == \n")
			
				for iZ = param.globalparam.N_iZ_Plot_Start:param.globalparam.N_iZ_Plot_End
					# << PLOT 1 >>
						Title = " iZ= $(IdSelect[iZ])"
						# Plot_∑infilt_Obs

							Label ="Obs_$(string(option.infilt.DataSingleDoubleRing⍰))_Ring"
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
							
						Path₂ = Path * "INFIL_" * string(option.infilt.Model⍰)  *  "_" * string(IdSelect[iZ]) *  ".svg"

					Plots.savefig(Plot_∑infilt_Obs, Path₂)
					println("    ~  $(Path₂) ~")
				end # for iZ=1:NiZ
			println("  ==  END: PLOT_∑INFILT  == \n")
			return nothing
			end # PLOT_∑INFILT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_∑INFILT_θΨ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PLOT_∑INFILT_θΨ(hydroInfilt, IdSelect, NiZ, optim, option, param, Path; hydro=[], N_Se=100)
			println("  ==  START: PLOT_∑INFILT_θΨ  ==")

				θ_Infilt      = fill(0.0::Float64, (N_Se))
				θ_Obs         = fill(0.0::Float64, (N_Se))
				Kunsat_Infilt = fill(0.0::Float64, (N_Se))
				Kunsat_Obs    = fill(0.0::Float64, (N_Se))

				for iZ = param.globalparam.N_iZ_Plot_Start: param.globalparam.N_iZ_Plot_End	
					Ψ_θΨobs_Min = 10.0 ^ -2 # [mm]

					Ψ_θΨobs_Max = 200000.0 * 10.0 # [mm]

					Ψ = 10.0 .^ range(log(Ψ_θΨobs_Min), stop=log(Ψ_θΨobs_Max), length=N_Se)

					θ_θΨobs_Max = hydroInfilt.Φ[iZ] + 0.1

					if option.run.HydroLabθΨ⍰ ≠ "No" && "Ks" ∈ optim.ParamOpt
						K_Ψ_Max = max(hydroInfilt.Ks[iZ], hydro.Ks[iZ]) * 1.1
					else
						K_Ψ_Max = hydroInfilt.Ks[iZ] * 1.1
					end #  "Ks" ∈ optim.ParamOpt

					for iΨ = 1:N_Se
						θ_Infilt[iΨ] = wrc.Ψ_2_θDual(option.infilt,Ψ[iΨ], iZ, hydroInfilt)

						Kunsat_Infilt[iΨ] = kunsat.Ψ_2_KUNSAT(option.infilt, Ψ[iΨ], iZ, hydroInfilt)

						if option.run.HydroLabθΨ⍰ ≠ "No"
							θ_Obs[iΨ] = wrc.Ψ_2_θDual(option.infilt,Ψ[iΨ], iZ, hydro)

							if option.run.HydroLabθΨ⍰ ≠ "No" && "Ks" ∈ optim.ParamOpt
								Kunsat_Obs[iΨ] = kunsat.Ψ_2_KUNSAT(option.infilt, Ψ[iΨ], iZ, hydro)
							end # "Ks" ∈ optim.ParamOpt		
						end # option.run.HydroLabθΨ⍰ ≠ :No
					end # iΨ 

					#PLOT 1:  Plot_θ_Ψ
						# Plot_θ_Ψ: Simulated Infiltration
							X = Ψ[1:N_Se] .* cst.Mm_2_Cm
							Y = θ_Infilt[1:N_Se]
							Label = "Infilt"
							Plot_θ_Ψ = Plots.plot(X, Y, seriestype=:line, label=Label, color= :blue, lw=2)

						# Plot_θ_Ψ: Observed
						if option.run.HydroLabθΨ⍰ ≠ "No"
							X = Ψ[1:N_Se] .* cst.Mm_2_Cm
							Y = θ_Obs[1:N_Se]
							Label = "Obs"
							Plot_θ_Ψ = Plots.plot!(X ,Y, seriestype=:line, label=Label, color= :red, lw=2)
						end # option.run.HydroLabθΨ⍰ ≠ :No

						# Plot_θ_Ψ: General attributes
							Plots.xlabel!("\\psi [cm]")
							Plots.ylabel!(L"\theta \ [cm^3 cm^{-3}]")
							Plots.plot!(xlims =(10.0*Ψ_θΨobs_Min*cst.Mm_2_Cm, Ψ_θΨobs_Max*cst.Mm_2_Cm), ylims =(0.0, θ_θΨobs_Max), xscale= :log10, size=(800,400), legend=:bottomleft)

						# PLOT2: Kunsat
							# Plot_K_Ψ: Obs K_Ψ
							X = Ψ[1:N_Se] .* cst.Mm_2_Cm
							Y = Kunsat_Infilt[1:N_Se] .* cst.MmS_2_CmH
							Label = "Infilt"
							Plot_K_Ψ = Plots.plot(X, Y, seriestype=:line, label=Label, color= :blue, lw=2)

							# Plot_K_Ψ: Sim K_Ψ
							if option.run.HydroLabθΨ⍰ ≠ "No" && "Ks" ∈ optim.ParamOpt
								X = Ψ[1:N_Se] .* cst.Mm_2_Cm
								Y = Kunsat_Obs[1:N_Se] .* cst.MmS_2_CmH
								Label = "Obs"
								Plot_K_Ψ = Plots.plot!(X, Y, seriestype=:line, label=Label, color= :red, lw=2)
							end # "Ks" ∈ optim.ParamOpt

							# General attributes
								Plots.xlabel!("\\psi [cm]")
								Plots.ylabel!(L" K (\psi) \ [cm \ h^{-1}]")
								Plot_K_Ψ = Plots.plot!(xlims = (Ψ_θΨobs_Min*cst.Mm_2_Cm, Ψ_θΨobs_Max*cst.Mm_2_Cm), ylims = (10^-2.0, K_Ψ_Max * cst.MmS_2_CmH), xscale= :log10,  yscale= :log10, legend=:bottomleft, size=(800,400))

						Path₂ = Path * "Infilt_ThetaH_" * string(option.hydro.HydroModel⍰) * "_" *string(IdSelect[iZ]) * ".svg"     
						Plot = Plots.plot(Plot_θ_Ψ, Plot_K_Ψ)
						Plots.savefig(Plot, Path₂)
						println("    ~  $(Path₂) ~")
				end # iZ

			println("  ==  END: PLOT_∑INFILT_θΨ  == \n")
			return nothing
			end  # function: PLOT_∑INFILT_θΨ

	end # module infilt
	# ............................................................

end  # module plot

