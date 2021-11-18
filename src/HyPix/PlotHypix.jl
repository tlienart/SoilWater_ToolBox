# =============================================================
#		module: plotHypix
# =============================================================
module plotHypix
	import  ..cst, ..kunsat, ..rootWaterUptake, ..tool, ..wrc, ..ΨminΨmax
	import Dates: value, DateTime
	# using PGFPlots

	# export θΨK

	# ========================================
	# PLOTTING HYDRAULIC RELATIONSHIP FOR EVERY HORIZON
	# ======================================== 
	# function θΨK(hydroHorizon, N_Layer, iOpt, pathHyPix)

	# 	# Deriving the Min and Max Ψ from principals of soil physics
	# 	Ψ_Min_Horizon = fill(0.0::Float64, N_Layer)
	# 	Ψ_Max_Horizon = fill(0.0::Float64, N_Layer)
	# 	for iZ=1:N_Layer
	# 		Ψ_Max_Horizon[iZ], Ψ_Min_Horizon[iZ] = ΨminΨmax.ΨMINΨMAX(hydroHorizon.θs[iZ], hydroHorizon.θsMacMat[iZ], hydroHorizon.σ[iZ], hydroHorizon.σMac[iZ], hydroHorizon.Ψm[iZ], hydroHorizon.ΨmMac[iZ])
	# 	end  # for iZ=1:N_Layer
		
	# 	# PREPARING THE DATA
	# 		N_Se = 1000
	# 		local Ψplot = exp.(range(log(minimum(Ψ_Min_Horizon[1:N_Layer])), stop = log(maximum(Ψ_Max_Horizon[1:N_Layer])), length=N_Se)) 

	# 		local θplot    = fill(0.0::Float64, N_Se)
	# 		local Kplot    = fill(0.0::Float64, N_Se)
	# 		local ∂θ∂Ψplot = fill(0.0::Float64, N_Se)
	# 		local ∂K∂Ψplot = fill(0.0::Float64, N_Se)

	# 		Plot_θΨK = PGFPlots.GroupPlot(4, 100, groupStyle = "horizontal sep = 3.5cm, vertical sep = 3.5cm")

	# 	# FOR EVERY HORIZON
	# 	for iZ = 1:N_Layer
			
	# 		for iΨ = 1:N_Se
	# 			if Ψ_Max_Horizon[iZ] ≥ Ψplot[iΨ] ≥ Ψ_Min_Horizon[iZ]
	# 				θplot[iΨ]    = wrc.Ψ_2_θDual(optionₘ,Ψplot[iΨ], iZ, hydroHorizon)
					
	# 				Kplot[iΨ]    = kunsat.Ψ_2_KUNSAT(optionₘ, Ψplot[iΨ], iZ, hydroHorizon)
					
	# 				∂θ∂Ψplot[iΨ] = wrc.∂θ∂Ψ(optionₘ, Ψplot[iΨ], iZ, hydroHorizon)

	# 				∂K∂Ψplot[iΨ] = kunsat.∂K∂Ψ(optionₘ, Ψplot[iΨ], iZ, hydroHorizon)
	# 			else
	# 				θplot[iΨ]    = NaN
					
	# 				Kplot[iΨ]    = NaN
					
	# 				∂θ∂Ψplot[iΨ] = NaN

	# 				∂K∂Ψplot[iΨ] = NaN
	# 			end
	# 		end # for iΨ

	# 		Θs_Max = maximum(hydroHorizon.θs[1:N_Layer]) + 0.05
	# 		Ks_Min = 10.0 ^ -7 * cst.MmS_2_CmH
	# 		Ks_Max = maximum(hydroHorizon.Ks[1:N_Layer]) * cst.MmS_2_CmH * 1.1

	# 		Title =" $(pathHyPix.SiteName_Hypix)  Layer = $(iZ)"
		
	# 	# Plot 1: θΨ
	# 		Plot_θΨ = PGFPlots.Plots.Linear(log.(Ψplot) , θplot, style=" smooth, blue, very thick", mark="none", legendentry=L"$ \theta ( \Psi ) $")

	# 		Plot_hydro = [Plot_θΨ]

	# 		push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xlabel=L"$ Ln \ \Psi [mm]$", ylabel=L"$ \theta \ [mm{^3} \ mm^{-3}]$", xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), ymin=0.0, ymax=Θs_Max, legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

	# 	# Plot 2: Kplot(Ψplot)
	# 		Plot_Kθ = PGFPlots.Plots.Linear(log.(Ψplot), Kplot .* cst.MmS_2_CmH, style=" smooth, red, very thick", mark="none", legendentry=L"$ K_{unsat} \ ( \Psi ) $")

	# 		Plot_hydro = [Plot_Kθ]
			
	# 		push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title,  xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), ymax=Ks_Max, ymode="log", xlabel=L"$Ln \  \Psi [mm]$", ylabel=L"$ K_{unsat} \ [cm \ h^{-1}]$", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

	# 	# Plot 3: ∂θ∂Ψplot
	# 		Plot_∂θ∂Ψ = PGFPlots.Plots.Linear(log.(Ψplot), ∂θ∂Ψplot , style=" smooth, green, very thick", mark="none", legendentry=L"$ \partial \theta \partial \Psi $")

	# 		Plot_hydro = [Plot_∂θ∂Ψ]

	# 		push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), xlabel=L"$Ln \  \Psi [mm] $", ylabel=L"$ \partial \theta \partial \Psi $", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

	# 	# Plot 4: ∂K∂Ψplot
	# 		Plot_∂K∂Ψ = PGFPlots.Plots.Linear(log.(Ψplot), ∂K∂Ψplot, style=" smooth, teal, very thick", mark="none", legendentry=L"$ \partial K \partial \Psi $")

	# 		Plot_hydro = [Plot_∂K∂Ψ]

	# 		push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), xlabel=L"$Ln \  \Psi \ [mm]$", ylabel=L"$ \partial K \partial \Psi $", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

	# 	end #iZ ............................................

	# 	Path = pathHyPix.Plot_Hypix_θΨK * "_" * string(iOpt) * ".svg"
	# 	PGFPlots.save(Path, Plot_θΨK) 
	# end # function θΨK


	# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# #		FUNCTION : ROOTDENSITY
	# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# 	function VEG_FUNCTIONS(discret, iOpt, N_iRoot, veg, Z, ΔRootDensity, pathHyPix)

	# 		Plot_All = PGFPlots.GroupPlot(2, 1, groupStyle = "horizontal sep = 3cm, vertical sep = 3cm")

	# 		# PLOT VEG_FUNCTIONS
	# 			ΔRootDensity_Norm = fill(0.0::Float64, N_iRoot)
	# 			# Taking accoung the tickness of the discretisation
	# 			# for iZ=1:N_iRoot
	# 			# 		ΔRootDensity_Norm[iZ] = Z[N_iRoot] * ΔRootDensity[iZ] / discret.ΔZ[iZ]
	# 			# end

	# 			# Plotting
	# 				Plot_RootDensity = PGFPlots.Plots.Linear(ΔRootDensity[1:N_iRoot], discret.Znode[1:N_iRoot], style=" smooth, cyan, very thick", mark="none")

	# 				Plot = [Plot_RootDensity]

	# 				push!(Plot_All, PGFPlots.Axis(Plot, style="width=12cm, height=8cm", xlabel=L"$ \Delta Rdf \ [\%] $", ylabel=L"$Z \ [mm]$", title="(a)"))
			
	# 		# PLOT StressReduction
	# 			# Data	
	# 			N_Se = 6
	# 			Ψstress = fill(0.0::Float64, 2, N_Se) 
	# 			Ψstress[1,1] = veg.Ψfeddes1 / 10.0
	# 			Ψstress[1,2] = veg.Ψfeddes1
	# 			Ψstress[1,3] = veg.Ψfeddes2
	# 			Ψstress[1,4] = veg.Ψfeddes3
	# 			Ψstress[1,5] = veg.Ψfeddes4
	# 			Ψstress[1,6] = veg.Ψfeddes4 * 2.0

	# 			Wsf = fill(0.0::Float64, N_Se)
	# 			for iΨ ∈ 1:N_Se
	# 				Wsf[iΨ] = rootWaterUptake.stressReduction.WATER_STRESS_FUNCTION(2, iΨ, veg, Ψstress)
	# 			end

	# 			Plot_Wsf = PGFPlots.Plots.Linear(Ψstress[1,1:N_Se] .* cst.Mm_2_kPa, Wsf[1:N_Se], style="violet, very thick", mark="none")

	# 			Plot = [Plot_Wsf]

	# 			push!(Plot_All, PGFPlots.Axis(Plot, style="width=12cm, height=8cm", xlabel=L"$ \Psi \ [kPa]$", xmode="log", ylabel=L"$ F_{waterStress} \ [-]$", title="(b)"))

	# 		Path = pathHyPix.Vegetation * "_" * string(iOpt) * ".svg"
	# 		PGFPlots.save(Path, Plot_All)	
	# 	end  # function ROOTDENSITY

	# =============================================================
	#		module: name
	# =============================================================
		module makkie
			using CairoMakie, LaTeXStrings
			using Dates
			export θPROFILE, TIMESERIES

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : θPROFILE(∑T_Reduced, discret, obsTheta, option, param, θ_Reduced)
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θPROFILE(∑T_Reduced, discret, iSim, NiZ, obsTheta, option, param, pathHyPix, Soilname, θ_Reduced)

				# PATH
					Path = pathHyPix.Plot_θprofile * ".svg"
					println("			 ~ ", Path, "~")
					rm(Path, force=true, recursive=true)


				# DEPTHS TO PLOT
					Zprofile = fill(0.0::Float64, NiZ)
					for iZ=1:NiZ
						Zprofile[iZ] = discret.Znode[iZ] / 100.0
					end

				# SELECTING PROFILE TO PLOT
					Nt = obsTheta.Nit

				# INITIALIZING PLOT
					CairoMakie.activate!()
					Makie.inline!(true)

					Color_Hypix = [:red, :darkviolet, :orange,  :blue, :teal]

					Fig = Figure(resolution=(600,500))
					Title = Soilname[iSim]

					Label_HyPix =fill("", Nt)
					Label_Hydrus =fill("", Nt)

					Ax1 = Axis(Fig[1,1], title=Title, xlabel= L"$\theta$  $[m^{3}  m^{-3}]$", ylabel= L"Z  $[cm]$",  font="CMU Serif", titlesize=25, fontsize=16, xlabelsize=22, ylabelsize=22, xgridvisible=false, ygridvisible=false)

					# Ax2 = Axis(Fig[1,1],  font = "CMU Serif", titlesize=30, fontsize=16, xlabelsize=24, ylabelsize=24)

				# For every θprofile_Time
					for iT=1:Nt
						Tprofile = obsTheta.∑T[iT]

						iTprofile = 1
		
						iTprofile = findfirst(x->x==Tprofile, ∑T_Reduced)
						
						if isnothing(iTprofile)
							println("Error θprofile_Time must be one of = $(obsTheta.∑T)")
							error()
						end

						θprofile = θ_Reduced[iTprofile, 1:NiZ]

						# PLOTTING

							# Label
								if param.hyPix.ΔT_Output==3600.0 
									Label_HyPix[iT] = "HyP=" * string(ceil(Int, Tprofile / param.hyPix.ΔT_Output)) * "Hour" 
								else
									Label_HyPix[iT] = "HyP_" * string(ceil(Int, Tprofile / param.hyPix.ΔT_Output)) * "Day" 
								end

								if param.hyPix.ΔT_Output==3600.0 
									Label_Hydrus[iT] = "HYD_" * string(ceil(Int, Tprofile / param.hyPix.ΔT_Output)) * "Hour" 
								else
									Label_Hydrus[iT] = "HYD_" * string(ceil(Int, Tprofile / param.hyPix.ΔT_Output)) * "Day" 
								end
	
						Plot2 = lines!(Ax1, obsTheta.θobs[iT,1:NiZ], -Zprofile, color=Color_Hypix[iT], linewidth=3, label=Label_Hydrus[iT])
						Plot1 = lines!(Ax1, θprofile, -Zprofile, color=Color_Hypix[iT], linewidth=2, label=Label_HyPix[iT], linestyle=:dash)
						
					end

					Leg = Legend(Fig[2,1], Ax1, framevisible=true, orientation=:horizontal, tellheight=true, nbanks=2, framecolor = (:grey,0.5), labelsize=14)
			
					trim!(Fig.layout)

					# axislegend()
					display(Fig)
					save(Path, Fig)


			return nothing
			end  # function: θPROFILE
			# ------------------------------------------------------------------

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : TIMESERIES
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function TIMESERIES(Date_Reduced, ∑T_Reduced, obsTheta, discret, iOpt, Nit_Reduced, NiZ, option, param, ΔEvaporation_Reduced, ΔQ_Reduced, ΔPrGross_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔSink_Reduced, θ_Reduced, θobs_Reduced, clim, i∑T_CalibrStart_Day, θsim_Aver, pathHyPix)

				# DAYS PLOT
					ΔDays = 60
				# OPTIONS
					Option_θobs = false

				# PATH
					Path = pathHyPix.Plot_HypixTime * "_" * string(iOpt) * ".svg"
					rm(Path, force=true, recursive=true)

				# TICKS
					# Date_Start_Calibr = obsTheta.Date[1]
					Date_Start_Calibr = obsTheta.Date[1]  # since we need to compute the culmulativeof the 1rst day
					
					Date_End_Calibr = obsTheta.Date[end] 
					
					DateTick = Date_Start_Calibr:Day(61):Date_End_Calibr
					
					DateTick2= Dates.format.(DateTick, "d u Y")
				
				# PLOTTING
				# , resolution = (3000, 2500)
					Fig = Figure( font="Sans", titlesize=40, fontsize=16, xlabelsize=24, ylabelsize=18)
				# Plot Climate	
				iSubplot = 0
				if option.hyPix.Plot_Climate
					iSubplot += 1
					Axis1 = Axis(Fig[iSubplot,1], title=pathHyPix.IdName_Hypix,  ylabel= L"$\Delta Fluxes$ $[mm$ $day ^{-1}]$", rightspinevisible = false)

					hidexdecorations!(Axis1, ticks=false, grid=false)
						
					Label1= L" $\Delta Pr$"
					Plot_Climate1 = barplot!(Axis1,  ∑T_Reduced[1:Nit_Reduced], ΔPrGross_Reduced[1:Nit_Reduced], strokecolor=:blue, strokewidth=1.5, color=:blue)
					
					Label2=L"$\Delta PrThrough$"					
					Plot_Climate2 = barplot!(Axis1, ∑T_Reduced[1:Nit_Reduced], ΔPr_Reduced[1:Nit_Reduced], strokecolor=:cyan, strokewidth=1, color=:cyan)
					
					Label3=L"$\Delta Hpond$"
					Plot_Climate3 = barplot!(Axis1, ∑T_Reduced[1:Nit_Reduced], ΔPond_Reduced[1:Nit_Reduced], strokecolor=:grey, strokewidth=1, colour=:grey)

					# ===
					Axis2  = Axis(Fig[iSubplot, 1], yticklabelcolor=:black, yaxisposition = :right, rightspinecolor = :black, ytickcolor=:black, ylabel= L"$\Delta EvapoTransp$ $[mm$ $day ^{-1}]$")

					hidexdecorations!(Axis2, ticks=false, grid=false)

					Label4 =L"$\Delta Pet$"
					Plot_Climate4 = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], ΔPet_Reduced[1:Nit_Reduced], linewidth=2, colour=:darkgreen)

					Label5 = L"$\Delta Sink$"
					Plot_Climate5 = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], ΔSink_Reduced[1:Nit_Reduced], linewidth=2, colour=:red)

					Label6=L"$\Delta Evap$"
					Plot_Climate6 = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], ΔEvaporation_Reduced[1:Nit_Reduced], linewidth=2, colour=:purple4)

					# Plot_Climate = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], (ΔSink_Reduced[1:Nit_Reduced].-ΔEvaporation_Reduced[1:Nit_Reduced]), colour=:blue, label=L"$\Delta Rwu$")
					
					# ==
										
					trim!(Fig.layout)

					iSubplot += 1
					Axis3  = Axis(Fig[iSubplot, 1], ylabel= L"$\Delta Q$ $[mm$ $day ^{-1}]$")
					hidexdecorations!(Axis3, ticks=false, grid=false)
					Plot_Climate7 = barplot!(Axis3, ∑T_Reduced[1:Nit_Reduced], -ΔQ_Reduced[1:Nit_Reduced, NiZ+1], strokecolor=:red, strokewidth=1, color=:red)

					iSubplot += 1

					Legend(Fig[1:2,2], [Plot_Climate1, Plot_Climate2, Plot_Climate3, Plot_Climate4, Plot_Climate5, Plot_Climate6], [Label1, Label2, Label3, Label4, Label5, Label6], framevisible=true, orientation=:vertical, tellwidth = false, tellheight = true)

					iSubplot += 1

					trim!(Fig.layout)
				end # if: option.hyPix.Plot_Climate

				# PLOT Θ
				if option.hyPix.Plot_θ
					
					Style_Hypix = [:red, :darkviolet, :orange, :teal, :blue]
					
					iSubplot += 1
					Axis4 = Axis(Fig[iSubplot,1], title=pathHyPix.IdName_Hypix, ylabel=L"$\theta$ $[mm^3 mm^{-3}]$")

					# Putting a monthly dates
						Ndates = length(∑T_Reduced)

						Ndates_Reduced = floor(Int, ∑T_Reduced[Ndates] / (ΔDays * param.hyPix.ΔT_Output)) + 1
						Date_Reduced2 = fill(Date_Reduced[1], Ndates_Reduced)
						Date_Reduced3 = fill("", Ndates_Reduced)
						∑T_Reduced2 = fill(0, Ndates_Reduced)

					iGood = 1
					for i=1:length(∑T_Reduced)

						if ∑T_Reduced[i] % (ΔDays * param.hyPix.ΔT_Output) == 0
							∑T_Reduced2[iGood] = ∑T_Reduced[i]	
							Date_Reduced2[iGood] = Date_Reduced[i]
							Date_Reduced3[iGood]= Dates.format(Date_Reduced2[iGood], "d u Y")
							iGood += 1
						end
					end

					# DateNewFormat = string.(Dates.format.(Date_Reduced2, "d u Y"))
					Axis4.xticks = (∑T_Reduced2[1:iGood],string.(Date_Reduced3[1:iGood]))
					Axis4.xticklabelrotation = π/4

					# Observation θplot obs
					for iZobs = 1:obsTheta.Ndepth
						# lABEL
							Label_Obs = "Obs=" * string(Int(floor(obsTheta.Z[iZobs]))) * "mm"

							Label_Sim = "θ" * string(Int(floor((discret.Znode[obsTheta.ithetaObs[iZobs]])))) * "mm"

							if Option_θobs
								Plot_θ = lines!(Axis4, ∑T_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, iZobs], linewidth=5, colour=Style_Hypix[iZobs], label=Label_Obs)
							end


							Plot_θ = lines!(Axis4, ∑T_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced, obsTheta.ithetaObs[iZobs]], linewidth=3, colour=Style_Hypix[iZobs], label=Label_Sim)
					end # loop

					iSubplot += 1
					Fig[iSubplot, 1] = Legend(Fig, Axis4, framevisible=true, orientation=:horizontal, tellheight=true, tellwidth = true, haligns=:center, valigns=:bottom)

					colgap!(Fig.layout, 90)
					rowgap!(Fig.layout, 5)
					trim!(Fig.layout)


						# Plot = plot(Plot, Plot_θ, Plot_Climate, xmin=Date_Reduced[1], xmax=Date_Reduced[Nit_Reduced], ymin=0.0, xtick=(DateTick,DateTick2), xrotation=rad2deg(pi/4), framestyle=:box, grid=true)

				end # if: option.hyPix.Plot_θ
				display(Fig)
				save(Path, Fig)
				println("			 ~ ", Path, "~")
			
			return nothing
			end  # function: TIMESERIES
			
		end  # module: makkie
		# ............................................................



			# =============================================================
			#		module: plots
			# =============================================================
			# module plots
			# import ...sorptivity, ..wrc, ..cst, ...reading
			# export PLOT_SORPTIVITY

			# 	using Plots.PlotMeasures, LaTeXStrings
			# 	using Plots
			# 	using Dates
				
			# 	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# 	#		FUNCTION : PLOT_SORPTIVITY
			# 	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# 	function PLOT_SORPTIVITY(hydro, iOpt, option, pathHyPix)
			# 		println("  ==  START: PLOT_SORPTIVITY_SeIni  ==")

			# 		# Setting the range of values for Se
         #          Se_Ini         = collect(0.0:0.001:1.0)
         #          N_SeIni        = length(Se_Ini)
         #          Sorptivity_Mod = fill(0.0::Float64, (N_SeIni))
         #          θini          = fill(0.0::Float64, (N_SeIni))

			# 		for iSeIni=1:N_SeIni
			# 			θini[iSeIni] = wrc.Se_2_θ(Se_Ini[iSeIni], 1, hydro)

			# 			Sorptivity_Mod[iSeIni] = sorptivity.SORPTIVITY(θini[iSeIni], 1, hydro, option) 
			# 		end
					
			# 		# PLOTTING ====================	
			# 			Plot1=Plots.plot(layout=1)

			# 			Title =" $(pathHyPix.SiteName_Hypix)"

			# 			Plots.plot!(Plot1, Se_Ini[1:N_SeIni] , Sorptivity_Mod[1:N_SeIni], framestyle = [:box :semi :origin :zerolines :grid :true], xlabel=L"Initial \ Se \ [-]", ylabel=L"Sorptivity \  [ \ mm \ \sqrt s \ ]", label="", grid=false) 
					
			# 			Path =pathHyPix.Plot_Sorptivity  * "_" * string(iOpt) * ".svg"

			# 			Plots.savefig(Plot1, Path)

			# 			println("			 ~ ", Path, "~")

			# 	end  # function: PLOT_SORPTIVITY


			# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# #		FUNCTION : TIMESERIES
			# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# 	function TIMESERIES(Date_Reduced, ∑T_Reduced, obsTheta, discret, iOpt, Nit_Reduced, NiZ, option, param, ΔEvaporation_Reduced, ΔQ_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔSink_Reduced, θ_Reduced, θobs_Reduced, clim, i∑T_CalibrStart_Day, θsim_Aver, pathHyPix)

			# 	# PATH
			# 		Path = pathHyPix.Plot_HypixTime * "_" * string(iOpt) * ".svg"
			# 		rm(Path, force=true, recursive=true)

			# 	# TICKS
			# 		# Date_Start_Calibr = obsTheta.Date[1]
			# 		Date_Start_Calibr = DateTime(param.hyPix.obsTheta.Year_Start, param.hyPix.obsTheta.Month_Start, param.hyPix.obsTheta.Day_Start, param.hyPix.obsTheta.Hour_Start, param.hyPix.obsTheta.Minute_Start, param.hyPix.obsTheta.Second_Start) # since we need to compute the culmulativeof the 1rst day
					
			# 		# Date_End_Calibr = obsTheta.Date[end]
			# 		Date_End_Calibr = DateTime(param.hyPix.Year_End, param.hyPix.Month_End, param.hyPix.Day_End, param.hyPix.Hour_End, param.hyPix.Minute_End, param.hyPix.Second_End)
					
			# 		DateTick=range(Date_Start_Calibr,step=Day(61),Date_End_Calibr)
					
			# 		DateTick2= Dates.format.(DateTick, "d u Y")
				
			# 	# PLOTTING
			# 		Plot = Plots.plot(layout=(3, 1), size=(2500,2200), bottom_margin=0.01mm)
					
			# 		default(titlefont=(20,"times"), legendfontsize=24, guidefont=18, tickfont=18, grid=true)

			# 	# Plot Climate	
			# 	iSubplot = 0
			# 	if option.hyPix.Plot_Climate
			# 		iSubplot += 1		

			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], -ΔQ_Reduced[1:Nit_Reduced, NiZ+1], label=L"$\Delta Q$", line=(:solid, 1), linecolour=:red, fillcolor=:darkred, fill=(0,:darkred))
					
			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], ΔPond_Reduced[1:Nit_Reduced], label=L"$\Delta H_{Pond}$", linecolour=:grey, fill = (0, :grey))
					
			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr[i∑T_CalibrStart_Day:clim.N_Climate], color=:blue, colorbar=false,  line =(:sticks, :solid, 5), label= L"$\Delta Pr  $")

			# 		Plot_Climate = Plots.plot!(Plot, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr_Through[i∑T_CalibrStart_Day:clim.N_Climate], color=:cyan, line =(:sticks, :solid, 4), colorbar=false, label=L"$\Delta Pr_through$")
	
			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ylabel=L"$Daily \ Simulation \ [mm]$", title=pathHyPix.IdName_Hypix, xtickfont = (0.01, :white), xrotation=rad2deg(pi/2))
			# 	end # if: option.hyPix.Plot_Climate

			# 	# PLOT EVAPOYTRANSPIRATION
			# 		iSubplot += 1	

			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[2:Nit_Reduced], ΔPet_Reduced[2:Nit_Reduced], linecolour=:darkgreen, label=L"$\Delta Pet$", line=(2.5,:solid))

			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[2:Nit_Reduced], ΔSink_Reduced[2:Nit_Reduced], linecolour=:red, line=(2.0,:solid), label=L"$\Delta Sink$")

			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[2:Nit_Reduced], (ΔSink_Reduced[2:Nit_Reduced].-ΔEvaporation_Reduced[2:Nit_Reduced]), label=L"$\Delta Rwu$", linecolour=:blue, line=(2.0,:solid))
					
			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[2:Nit_Reduced], ΔEvaporation_Reduced[2:Nit_Reduced], label=L"$\Delta Evap$", linecolour=:purple4, line=(2.0,:solid))

			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ylabel=L"$Daily \ Simulation \ [mm]$", xtickfont = (0.01, :white), xrotation=rad2deg(pi/2))

			# 	# PLOT Θ
			# 	if option.hyPix.Plot_θ
			# 		iSubplot += 1

			# 		Style_Hypix = [:red, :darkviolet, :orange, :teal, :blue]

			# 		# Observation θplot obs
			# 		for ithetaObs = 1:obsTheta.Ndepth
			# 			# lABEL
			# 				Label_Obs = "Obs=" * string(Int(floor(obsTheta.Z[ithetaObs]))) * "mm"

			# 				Label_Sim = "Sim=" * string( Int(floor((discret.Znode[obsTheta.ithetaObs[ithetaObs]])))) * "mm"

			# 			# Plotting
			# 				# Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, iZobs].+param.hypixStart.calibr.θobs_Uncert, line=(0.5,:solid), linecolour=Style_Hypix[iZobs], label=false)
		
			# 				# Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], max.(θobs_Reduced[1:Nit_Reduced, iZobs].-param.hypixStart.calibr.θobs_Uncert, 0.0), line=(0.5,:solid), linecolour=Style_Hypix[iZobs], label=false)
							
			# 				if option.hyPix.θobs_Average
			# 					Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, ithetaObs], line=(2.5,:solid), linecolour=Style_Hypix[ithetaObs], label="Obs θaver [0-40cm]")

			# 					Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θsim_Aver[1:Nit_Reduced], label="Sim θaver [0-40cm]", line=(2.5,:solid), linecolour=:blue)

			# 					Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced,4], label="Sim θ=10cm", line=(2.5,:dashdot), linecolour=:darkblue)

			# 					Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced,14], label="Sim θ=35cm", line=(2.5,:dashdot), linecolour=:darkblue)
			# 			else
			# 				Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, iZobs], line=(2.5,:solid), linecolour=Style_Hypix[iZobs], label=Label_Obs)

			# 				Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced, calibr.iZobs[iZobs]], label=Label_Sim, line=(2.5,:dashdot), linecolour=Style_Hypix[iZobs])
			# 			end  # if: option.hyPix.

			# 		end # loop

			# 		Plot_θ = Plots.plot!(subplot=iSubplot, ylabel=L"$\theta \ [mm^3 \ mm^{-3}]$")

			# 		Plot = Plots.plot(Plot, Plot_θ, Plot_Climate, xmin=Date_Reduced[1], xmax=Date_Reduced[Nit_Reduced], ymin=0.0, xtick=(DateTick,DateTick2), xrotation=rad2deg(pi/4), framestyle=:box, grid=true)

			# 	end # if: option.hyPix.Plot_θ
				
			# 	Plots.savefig(Plot, Path)
			# 	println("			 ~ ", Path, "~")
			
			# 	return nothing
			# 	end  # function: TIMESERIES

			# end  # module: plots
			# ............................................................

end  # module plotHypix
# ............................................................