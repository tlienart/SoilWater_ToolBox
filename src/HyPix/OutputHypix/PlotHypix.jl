# =============================================================
#		module: plotHypix
# =============================================================
module plotHypix
	import  ..cst, ..kunsat, ..rootWaterUptake, ..tool, ..wrc, ..ΨminΨmax
	import Dates: value, DateTime


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT_HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_HYPIX(∑T_Reduced, clim, Date_Reduced, discret, hydro, hydroHorizon, i∑T_CalibrStart_Day, iMultistep, iScenario, N_iRoot, N_Layer, Nit_Reduced, NiZ, obsTheta, optionHypix, paramHypix, pathOutputHypix, SiteName, veg, Z, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔPrGross_Reduced, ΔQ_Reduced, ΔRootDensity, ΔSink_Reduced, θ_Reduced, θobs_Reduced, θsim_Aver)

			if optionHypix.Plot_Hypix
				plotHypix.makkie.TIMESERIES(∑T_Reduced, clim, Date_Reduced, discret, i∑T_CalibrStart_Day, iMultistep, iScenario, Nit_Reduced, NiZ, obsTheta, optionHypix, paramHypix,  pathOutputHypix, SiteName, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔPrGross_Reduced, ΔQ_Reduced, ΔSink_Reduced, θ_Reduced, θobs_Reduced, θsim_Aver)
			end

			if optionHypix.Plot_θprofile
				plotHypix.makkie.θPROFILE(∑T_Reduced, discret, iScenario, NiZ, obsTheta, optionHypix, paramHypix, pathOutputHypix, SiteName, θ_Reduced)
			end  # if: optionHypix.Plot_
			if optionHypix.Plot_θΨK
				plotHypix.θΨK(hydroHorizon, N_Layer, iMultistep, pathOutputHypix)
			end
			if optionHypix.Plot_Vegetation && optionHypix.RootWaterUptake
				plotHypix.VEG_FUNCTIONS(discret, iMultistep, N_iRoot, veg, Z, ΔRootDensity, pathOutputHypix)
			end
			if optionHypix.Plot_Interception
				plotHypix.plots.RAINFALL_INTERCEPTION(clim, i∑T_CalibrStart_Day, iMultistep, pathOutputHypix)
			end
			if  optionHypix.Plot_Sorptivity
				plotHypix.plots.PLOT_SORPTIVITY(hydro, iMultistep, optionHypix, pathOutputHypix)
			end
			

		return nothing
		end  # function: PLOT_HYPIX
	# ------------------------------------------------------------------


	# export θΨK

	# ========================================
	# PLOTTING HYDRAULIC RELATIONSHIP FOR EVERY HORIZON
	# ======================================== 
	# function θΨK(hydroHorizon, N_Layer, iMultistep, pathHyPix)

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

	# 				∂K∂Ψplot[iΨ] = kunsat.∂K∂ΨMODEL(optionₘ, Ψplot[iΨ], iZ, hydroHorizon)
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

	# 	Path = pathHyPix.Plot_Hypix_θΨK * "_" * string(iMultistep) * ".svg"
	# 	PGFPlots.save(Path, Plot_θΨK) 
	# end # function θΨK




	# =============================================================
	#		module: makkie
	# =============================================================
		module makkie
			using CairoMakie
			using Dates
			export θPROFILE, TIMESERIES

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : θPROFILE(∑T_Reduced, discret, obsTheta, optionHypix, paramHypix, θ_Reduced)
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θPROFILE(∑T_Reduced, discret, iScenario, NiZ, obsTheta, optionHypix, paramHypix, pathHyPix, SiteName, θ_Reduced)

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
					Title = SiteName[iScenario]

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
								if paramHypix.ΔT_Output==3600.0 
									Label_HyPix[iT] = "HyP=" * string(ceil(Int, Tprofile / paramHypix.ΔT_Output)) * "Hour" 
								else
									Label_HyPix[iT] = "HyP_" * string(ceil(Int, Tprofile / paramHypix.ΔT_Output)) * "Day" 
								end

								if paramHypix.ΔT_Output==3600.0 
									Label_Hydrus[iT] = "HYD_" * string(ceil(Int, Tprofile / paramHypix.ΔT_Output)) * "Hour" 
								else
									Label_Hydrus[iT] = "HYD_" * string(ceil(Int, Tprofile / paramHypix.ΔT_Output)) * "Day" 
								end
	
						Plot2 = lines!(Ax1, obsTheta.θobs[iT,1:NiZ], -Zprofile, color=Color_Hypix[iT], linewidth=3, label=Label_Hydrus[iT])
						Plot1 = lines!(Ax1, θprofile, -Zprofile, color=Color_Hypix[iT], linewidth=2, label=Label_HyPix[iT], linestyle=:dash)
						
					end

					Leg = Legend(Fig[2,1], Ax1, framevisible=true, orientation=:horizontal, tellheight=true, nbanks=2, labelsize=14)
			
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
			function TIMESERIES(∑T_Reduced::Vector{Float64}, clim, Date_Reduced::Vector{Dates.DateTime}, discret, i∑T_CalibrStart_Day::Int64, iMultistep::Int64, iScenario::Int64, Nit_Reduced::Int64, NiZ::Int64, obsTheta, optionHypix, paramHypix, pathHyPix, SiteName::Vector{Any}, ΔEvaporation_Reduced::Vector{Float64}, ΔPet_Reduced::Vector{Float64}, ΔPond_Reduced::Vector{Float64}, ΔPr_Reduced::Vector{Float64}, ΔPrGross_Reduced::Vector{Float64}, ΔQ_Reduced::Matrix{Float64}, ΔSink_Reduced::Vector{Float64}, θ_Reduced::Matrix{Float64}, θobs_Reduced::Matrix{Float64}, θsim_Aver::Vector{Float64})

				# Number of days to print 
					Nticks = 12

				# PATH
					Path = pathHyPix.Plot_HypixTime * "_" * string(iMultistep) * ".svg"
					rm(Path, force=true, recursive=true)

				# STYLE colour
					Style_Color = [:red, :darkviolet, :orange, :teal, :blue, :green]

				# TICKS
					# Date_Start_Calibr = obsTheta.Date[1]
					Date_Start_Calibr = obsTheta.Date[1]  # since we need to compute the culmulativeof the 1rst day
					
					Date_End_Calibr = obsTheta.Date[end]
					
				# DAYS PLOT
					ΔDays = floor(Int, Dates.value((Date_End_Calibr - Date_Start_Calibr)) / ((Nticks - 1) * paramHypix.ΔT_Output * 1000))::Int64

					# Putting a monthly dates
						Ndates = length(∑T_Reduced)

						Ndates_Reduced = floor(Int, ∑T_Reduced[Ndates] / (ΔDays * paramHypix.ΔT_Output)) + 1
						
						Date_Reduced2 = fill("", Ndates_Reduced+1)
						∑T_Reduced2 = fill(0::Int64, Ndates_Reduced+1)
	
						# Reducing dates
						∑T_Reduced2[1] = ∑T_Reduced[1]	
		
						Date_Reduced2[1]= Dates.format(Date_Reduced[1], "d u Y")
						iGood = 2
						for i=1:Ndates
							if ∑T_Reduced[i] -  ∑T_Reduced[1] ≥ ((iGood - 1) * ΔDays * paramHypix.ΔT_Output)

								∑T_Reduced2[iGood] = ∑T_Reduced[i]	
		
								Date_Reduced2[iGood]= Dates.format(Date_Reduced[i], "d u Y")

								iGood += 1
							end # if
						end # for
				
				# PLOTTING
				# , resolution = (3000, 2500)
					Fig = Figure( font="Sans", titlesize=40, fontsize=12, xlabelsize=13, ylabelsize=13, labelsize=12, resolution = (1000, 700))
				
				# Plot Climate -------------------------------------------------	
				iSubplot = 0
				if optionHypix.Plot_Climate
					iSubplot += 1

					Axis1 = Axis(Fig[iSubplot,1], title=SiteName[iScenario], ylabel= L"$\Delta Fluxes$ $[mm$ $day ^{-1}]$", rightspinevisible = false, xgridvisible=false, ygridvisible=false)

					xlims!(Axis1, ∑T_Reduced[1], ∑T_Reduced[Nit_Reduced])
					ylims!(Axis1, low=0)

					hidexdecorations!(Axis1, ticks=false)
						
					Label1= L" $\Delta Pr$"
					Plot_Climate1 = barplot!(Axis1, ∑T_Reduced[1:Nit_Reduced], ΔPrGross_Reduced[1:Nit_Reduced], strokecolor=:blue, strokewidth=1.5, color=:blue)
					
					Label2=L"$\Delta PrThrough$"					
					Plot_Climate2 = barplot!(Axis1, ∑T_Reduced[1:Nit_Reduced], ΔPr_Reduced[1:Nit_Reduced], strokecolor=:cyan, strokewidth=1, color=:cyan)
					
					Label3=L"$\Delta Hpond$"
					Plot_Climate3 = barplot!(Axis1, ∑T_Reduced[1:Nit_Reduced], ΔPond_Reduced[1:Nit_Reduced], strokecolor=:grey, strokewidth=1, colour=:grey)

					# ---------------------------------------------------------
					Axis2  = Axis(Fig[iSubplot, 1], yticklabelcolor=:black, yaxisposition = :right, rightspinecolor = :black, ytickcolor=:black, ylabel= L"$\Delta ET$ $[mm$ $day ^{-1}]$", xgridvisible=false, ygridvisible=false)
					xlims!(Axis2, ∑T_Reduced[1],∑T_Reduced[Nit_Reduced])
					ylims!(Axis2, low=0)

					hidexdecorations!(Axis2, ticks=false)

					Label4 =L"$\Delta Pet$"
					Plot_Climate4 = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], ΔPet_Reduced[1:Nit_Reduced], linewidth=2, colour=:darkgreen)

					Label5 = L"$\Delta Sink$"
					Plot_Climate5 = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], ΔSink_Reduced[1:Nit_Reduced], linewidth=2, colour=:red)

					Label6=L"$\Delta Evap$"
					Plot_Climate6 = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], ΔEvaporation_Reduced[1:Nit_Reduced], linewidth=2, colour=:purple4)

					# Plot_Climate = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], (ΔSink_Reduced[1:Nit_Reduced].-ΔEvaporation_Reduced[1:Nit_Reduced]), colour=:blue, label=L"$\Delta Rwu$")
					
					# GROUND WATER RECHARGE AXIS 3 ===================================================================										
						trim!(Fig.layout)

						Label7= L"$\Delta Q$"
						iSubplot += 1
						Axis3  = Axis(Fig[iSubplot, 1], ylabel= L"$\Delta Q$ $[mm$ $day ^{-1}]$", xgridvisible=false, ygridvisible=false)
						xlims!(Axis3, ∑T_Reduced[1],∑T_Reduced[Nit_Reduced])
						ylims!(Axis3, -maximum(ΔQ_Reduced[1:Nit_Reduced, NiZ+1]), 0)

						hidexdecorations!(Axis3)
						Plot_Climate7 = barplot!(Axis3, ∑T_Reduced[1:Nit_Reduced], -ΔQ_Reduced[1:Nit_Reduced, NiZ+1], strokecolor=:red, strokewidth=1, color=:red)

						iSubplot += 1

						Legend(Fig[iSubplot,1], [Plot_Climate1, Plot_Climate2, Plot_Climate3, Plot_Climate7, Plot_Climate4, Plot_Climate5, Plot_Climate6], [Label1, Label2, Label3, Label7, Label4, Label5, Label6], framevisible=true, orientation=:horizontal, tellwidth=true, nbanks=1)

						iSubplot += 1

						trim!(Fig.layout)
				end # if: optionHypix.Plot_Climate

				# PLOT Θ
					if optionHypix.Plot_θ		
						iSubplot += 1
						Axis4 = Axis(Fig[iSubplot,1], ylabel=L"$\theta$ $[mm^3 mm^{-3}]$", xgridvisible=false, ygridvisible=false)
						xlims!(Axis4, ∑T_Reduced[1],∑T_Reduced[Nit_Reduced])
						
							Axis4.xticks = (∑T_Reduced2[1:iGood], string.(Date_Reduced2[1:iGood]))
							Axis4.xticklabelrotation = π/4

							# Observation θplot obs
							for iZobs = 1:obsTheta.Ndepth
								# lABEL
									Label_Obs = "θobs" * string(Int(floor(obsTheta.Z[iZobs]))) * "mm"

									Label_Sim = "θsim" * string(Int(floor((discret.Znode[obsTheta.ithetaObs[iZobs]])))) * "mm"

									if optionHypix.θobs
										Plot_θobs = lines!(Axis4, ∑T_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, iZobs], linewidth=1.5, color=Style_Color[iZobs], label=Label_Obs)
									end

									Plot_θsim = lines!(Axis4, ∑T_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced, obsTheta.ithetaObs[iZobs]], linewidth=1.5, color=Style_Color[iZobs], label=Label_Sim, linestyle = :dash)
							end # loop

							iSubplot += 1
							# Fig[iSubplot, 1] = Legend(Fig, Axis4, framevisible=true, orientation=:horizontal, tellwidth = true, haligns=:center, valigns=:bottom, nbanks=2)

						colgap!(Fig.layout, 90)
						rowgap!(Fig.layout, 10)
						trim!(Fig.layout)

				end # if: optionHypix.Plot_θ
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
			# import ...sorptivity, ..wrc, ..cst
			# export PLOT_SORPTIVITY

			# 	using Plots.PlotMeasures, LaTeXStrings
			# 	using Plots
			# 	using Dates
				
			# 	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# 	#		FUNCTION : PLOT_SORPTIVITY
			# 	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# 	function PLOT_SORPTIVITY(hydro, iMultistep, optionHypix, pathHyPix)
			# 		println("  ==  START: PLOT_SORPTIVITY_SeIni  ==")

			# 		# Setting the range of values for Se
         #          Se_Ini         = collect(0.0:0.001:1.0)
         #          N_SeIni        = length(Se_Ini)
         #          Sorptivity_Mod = fill(0.0::Float64, (N_SeIni))
         #          θini          = fill(0.0::Float64, (N_SeIni))

			# 		for iSeIni=1:N_SeIni
			# 			θini[iSeIni] = wrc.Se_2_θ(Se_Ini[iSeIni], 1, hydro)

			# 			Sorptivity_Mod[iSeIni] = sorptivity.SORPTIVITY(θini[iSeIni], 1, hydro, optionHypix) 
			# 		end
					
			# 		# PLOTTING ====================	
			# 			Plot1=Plots.plot(layout=1)

			# 			Title =" $(pathHyPix.SiteName_Hypix)"

			# 			Plots.plot!(Plot1, Se_Ini[1:N_SeIni] , Sorptivity_Mod[1:N_SeIni], framestyle = [:box :semi :origin :zerolines :grid :true], xlabel=L"Initial \ Se \ [-]", ylabel=L"Sorptivity \  [ \ mm \ \sqrt s \ ]", label="", grid=false) 
					
			# 			Path =pathHyPix.Plot_Sorptivity  * "_" * string(iMultistep) * ".svg"

			# 			Plots.savefig(Plot1, Path)

			# 			println("			 ~ ", Path, "~")

			# 	end  # function: PLOT_SORPTIVITY


			# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# #		FUNCTION : TIMESERIES
			# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# 	function TIMESERIES(Date_Reduced, ∑T_Reduced, obsTheta, discret, iMultistep, Nit_Reduced, NiZ, optionHypix, paramHypix, ΔEvaporation_Reduced, ΔQ_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔSink_Reduced, θ_Reduced, θobs_Reduced, clim, i∑T_CalibrStart_Day, θsim_Aver, pathHyPix)

			# 	# PATH
			# 		Path = pathHyPix.Plot_HypixTime * "_" * string(iMultistep) * ".svg"
			# 		rm(Path, force=true, recursive=true)

			# 	# TICKS
			# 		# Date_Start_Calibr = obsTheta.Date[1]
			# 		Date_Start_Calibr = DateTime(paramHypix.opt.Year_Start, paramHypix.opt.Month_Start, paramHypix.opt.Day_Start, paramHypix.opt.Hour_Start, paramHypix.opt.Minute_Start, paramHypix.opt.Second_Start) # since we need to compute the culmulativeof the 1rst day
					
			# 		# Date_End_Calibr = obsTheta.Date[end]
			# 		Date_End_Calibr = DateTime(paramHypix.Year_End, paramHypix.Month_End, paramHypix.Day_End, paramHypix.Hour_End, paramHypix.Minute_End, paramHypix.Second_End)
					
			# 		DateTick=range(Date_Start_Calibr,step=Day(61),Date_End_Calibr)
					
			# 		DateTick2= Dates.format.(DateTick, "d u Y")
				
			# 	# PLOTTING
			# 		Plot = Plots.plot(layout=(3, 1), size=(2500,2200), bottom_margin=0.01mm)
					
			# 		default(titlefont=(20,"times"), legendfontsize=24, guidefont=18, tickfont=18, grid=true)

			# 	# Plot Climate	
			# 	iSubplot = 0
			# 	if optionHypix.Plot_Climate
			# 		iSubplot += 1		

			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], -ΔQ_Reduced[1:Nit_Reduced, NiZ+1], label=L"$\Delta Q$", line=(:solid, 1), linecolour=:red, fillcolor=:darkred, fill=(0,:darkred))
					
			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], ΔPond_Reduced[1:Nit_Reduced], label=L"$\Delta H_{Pond}$", linecolour=:grey, fill = (0, :grey))
					
			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr[i∑T_CalibrStart_Day:clim.N_Climate], color=:blue, colorbar=false,  line =(:sticks, :solid, 5), label= L"$\Delta Pr  $")

			# 		Plot_Climate = Plots.plot!(Plot, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr_Through[i∑T_CalibrStart_Day:clim.N_Climate], color=:cyan, line =(:sticks, :solid, 4), colorbar=false, label=L"$\Delta Pr_through$")
	
			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ylabel=L"$Daily \ Simulation \ [mm]$", title=pathHyPix.IdName_Hypix, xtickfont = (0.01, :white), xrotation=rad2deg(pi/2))
			# 	end # if: optionHypix.Plot_Climate

			# 	# PLOT EVAPOYTRANSPIRATION
			# 		iSubplot += 1	

			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[2:Nit_Reduced], ΔPet_Reduced[2:Nit_Reduced], linecolour=:darkgreen, label=L"$\Delta Pet$", line=(2.5,:solid))

			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[2:Nit_Reduced], ΔSink_Reduced[2:Nit_Reduced], linecolour=:red, line=(2.0,:solid), label=L"$\Delta Sink$")

			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[2:Nit_Reduced], (ΔSink_Reduced[2:Nit_Reduced].-ΔEvaporation_Reduced[2:Nit_Reduced]), label=L"$\Delta Rwu$", linecolour=:blue, line=(2.0,:solid))
					
			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[2:Nit_Reduced], ΔEvaporation_Reduced[2:Nit_Reduced], label=L"$\Delta Evap$", linecolour=:purple4, line=(2.0,:solid))

			# 		Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ylabel=L"$Daily \ Simulation \ [mm]$", xtickfont = (0.01, :white), xrotation=rad2deg(pi/2))

			# 	# PLOT Θ
			# 	if optionHypix.Plot_θ
			# 		iSubplot += 1

			# 		Style_Color = [:red, :darkviolet, :orange, :teal, :blue]

			# 		# Observation θplot obs
			# 		for ithetaObs = 1:obsTheta.Ndepth
			# 			# lABEL
			# 				Label_Obs = "Obs=" * string(Int(floor(obsTheta.Z[ithetaObs]))) * "mm"

			# 				Label_Sim = "Sim=" * string( Int(floor((discret.Znode[obsTheta.ithetaObs[ithetaObs]])))) * "mm"

			# 			# Plotting
			# 				# Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, iZobs].+paramHypix.hypixStart.calibr.θobs_Uncert, line=(0.5,:solid), linecolour=Style_Color[iZobs], label=false)
		
			# 				# Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], max.(θobs_Reduced[1:Nit_Reduced, iZobs].-paramHypix.hypixStart.calibr.θobs_Uncert, 0.0), line=(0.5,:solid), linecolour=Style_Color[iZobs], label=false)
							
			# 				if optionHypix.θavr_RootZone
			# 					Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, ithetaObs], line=(2.5,:solid), linecolour=Style_Color[ithetaObs], label="Obs θaver [0-40cm]")

			# 					Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θsim_Aver[1:Nit_Reduced], label="Sim θaver [0-40cm]", line=(2.5,:solid), linecolour=:blue)

			# 					Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced,4], label="Sim θ=10cm", line=(2.5,:dashdot), linecolour=:darkblue)

			# 					Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced,14], label="Sim θ=35cm", line=(2.5,:dashdot), linecolour=:darkblue)
			# 			else
			# 				Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, iZobs], line=(2.5,:solid), linecolour=Style_Color[iZobs], label=Label_Obs)

			# 				Plot_θ = Plots.plot!(Plot, subplot=iSubplot, Date_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced, calibr.iZobs[iZobs]], label=Label_Sim, line=(2.5,:dashdot), linecolour=Style_Color[iZobs])
			# 			end  # if: optionHypix.

			# 		end # loop

			# 		Plot_θ = Plots.plot!(subplot=iSubplot, ylabel=L"$\theta \ [mm^3 \ mm^{-3}]$")

			# 		Plot = Plots.plot(Plot, Plot_θ, Plot_Climate, xmin=Date_Reduced[1], xmax=Date_Reduced[Nit_Reduced], ymin=0.0, xtick=(DateTick,DateTick2), xrotation=rad2deg(pi/4), framestyle=:box, grid=true)

			# 	end # if: optionHypix.Plot_θ
				
			# 	Plots.savefig(Plot, Path)
			# 	println("			 ~ ", Path, "~")
			
			# 	return nothing
			# 	end  # function: TIMESERIES

			# end  # module: plots
			# ............................................................

end  # module plotHypix
# ............................................................