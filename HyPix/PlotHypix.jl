# =============================================================
#		module: plotHypix
# =============================================================
module plotHypix

	import  ..cst, ..kunsat, ..option, ..param, ..path, ..rootwateruptake, ..tool, ..wrc, ..ΨminΨmax
	import Dates: value, DateTime
	using PGFPlots

	export TIME_SERIES, RAINFALL_INTERCEPTION, θΨK

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_SERIES
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_SERIES(∑T_Plot, ∑WaterBalance_η_Plot, calibr, discret, Flag_Plot_Pond, iSim, N_∑T_Plot, N_iZ, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, ΔT_Plot, θ_Plot, θobs_Plot, Ψ_Plot)
				
			# PLOTTING GENERAL	
				Plot_TimeSeries = PGFPlots.GroupPlot(1, 10, groupStyle = "vertical sep = 3.5cm")

				Style_Obs = ["smooth, red, solid, ultra thick", "smooth, violet, solid, ultra thick", "smooth, orange, solid, ultra thick", "smooth, teal, ultra thick", "smooth, blue, ultra thick"]

				Style_Hypix = ["smooth, red, densely dashed, ultra thick", "smooth, violet, densely dashed, ultra thick", "smooth, orange, densely dashed, ultra thick", "smooth, teal, densely dashed, ultra thick", "smooth, blue, densely dashed, ultra thick"]

				CellTop=[1,16]
				Style_HypixTop = ["smooth, blue, densely dashed, thick", "smooth, cyan, densely dashed, thick"]

			# PLOTTING GRAPHs
				#  Plotting Pr and Pet
				if option.hypix.Plot_Climate
					Plot_Pr = PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], ΔPr_Plot[1:N_∑T_Plot].*cst.Mm_2_Cm, style="ybar,fill=blue, blue", mark="none", legendentry=L"$ Pr_{Through} \ [cm \ day^{-1}] $")

					Plot_Pet = PGFPlots.Plots.Linear(∑T_Plot[2:N_∑T_Plot], ΔPet_Plot[2:N_∑T_Plot], mark="none", style=" smooth, red, very thick", legendentry= L"$ Pet \ [mm \ day^{-1}] $")

					Plot_RootWaterUptake = PGFPlots.Plots.Linear(∑T_Plot[2:N_∑T_Plot], ΔSink_Plot[2:N_∑T_Plot], style=" smooth, teal, very thick, densely dashed", mark="none", legendentry=L"$ RootWaterUptake \ [mm \ day^{-1}] $")

					if Flag_Plot_Pond
						Plot_Pond = PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], ΔPond_Plot[1:N_∑T_Plot], style="ycomb, violet", mark="none", legendentry=L"$ Pond \ [mm \ \Delta T ^{-1}] $")

						Plot_Climate = [Plot_Pr; Plot_Pet; Plot_RootWaterUptake; Plot_Pond]
					else
						Plot_Climate = [Plot_Pr; Plot_Pet; Plot_RootWaterUptake]
					end # Flag_Plot_Pond	

					push!(Plot_TimeSeries, PGFPlots.Axis(Plot_Climate, style="width=20cm, height=8cm", xlabel=L"$\ Time \ [Day]$", ylabel=L"$Pr \ [cm \ day^{-1}]  \ \ Pet \ [mm \ day^{-1}]$", title=path.SiteName_Hypix, xmax=∑T_Plot[N_∑T_Plot], ymin=0.0, legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=4}"))
				end # if: option.hypix.Plot_Climate
					
				# PLOT Θ
				if option.hypix.Plot_θ
					# Observation θplot obs
						Plot_θobs = map(1:calibr.Ndepth) do iZobs
							Label_Obs = "Obs=" * string((calibr.Z[iZobs])) * "mm"

						PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], θobs_Plot[1:N_∑T_Plot, iZobs], legendentry=Label_Obs, style=Style_Obs[iZobs], mark="none")		  
						end # loop

					# Simulation θplot
						Plot_θhypix = map(1:calibr.Ndepth) do iZobs
							Label_Sim = "HyPix=" * string( (discret.Znode[calibr.iZobs[iZobs]])) * "mm"

							PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], θ_Plot[1:N_∑T_Plot, calibr.iZobs[iZobs]], style=Style_Hypix[iZobs], legendentry=Label_Sim, mark="none")
						end # loop

					# Simulation θtop
						Plot_θhypixTop = map(1:length(CellTop)) do iTop
							Label_HypixTop = "HyPix=" * string(discret.Znode[CellTop[iTop]]) * "mm"

							PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], θ_Plot[1:N_∑T_Plot,CellTop[iTop]] ,style=Style_HypixTop[iTop], legendentry=Label_HypixTop, mark="none")
						end # loop

					Ploting_θ = [Plot_θobs; Plot_θhypix]

					push!(Plot_TimeSeries, PGFPlots.Axis(Ploting_θ, style="width=20cm, height=8cm", xlabel=L"$\ Time \ [Day]$", ylabel=L"$\theta \ [mm^3 \ mm^{-3}]$", xmax=∑T_Plot[N_∑T_Plot], legendStyle ="{at={(0.0,-0.3)}, anchor=south west, legend columns=5}"))
					#
				end # if: option.hypix.Plot_θ

				# PLOTTING GRAPH WATER BALANCE
				if option.hypix.Plot_WaterBalance
					Plot_Wb = PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], ∑WaterBalance_η_Plot[1:N_∑T_Plot], style="smooth, red, very thick", mark="none", legendentry=L"$∑WaterBalance_η$")

					# Plot_Line = PGFPlots.Plots.Linear(x->1, (0, ∑T_Plot[N_∑T_Plot]), xbins=10, style="blue, dashed")

					Ploting = [Plot_Wb]

					push!(Plot_TimeSeries, PGFPlots.Axis(Ploting, style="width=20cm, height=8cm", xlabel=L"$\ Time \ [Day]$", ylabel=L"$ ∑WaterBalance_η \ [mm]$", xmax=∑T_Plot[N_∑T_Plot], legendStyle ="{at={(0.0,-0.3)}, anchor=south west, legend columns=1}"))
				end # if: Plot_WaterBalance


				# PLOTTING ΔT
				if option.hypix.Plot_ΔT
					Plot_ΔT = PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], ΔT_Plot[1:N_∑T_Plot], style="smooth, cyan, very thick", mark="none", legendentry=L"$ \Delta T $")

					Ploting_ΔT = [Plot_ΔT]

					push!(Plot_TimeSeries, PGFPlots.Axis(Ploting_ΔT, style="width=20cm, height=8cm", xlabel=L"$\ Time \ [Day]$", ylabel=L"$ \Delta T \ [mm]$",xmax=∑T_Plot[N_∑T_Plot], ymode="log", legendStyle ="{at={(0.0,-0.3)}, anchor=south west, legend columns=1}"))
				end

				
				# PLOTTING GRAPH Ψplot
				if option.hypix.Plot_Ψ
						Plot_Ψhypix = map(1:calibr.Ndepth) do iZobs
							Label_Sim = "HyPix=" * string( (discret.Znode[calibr.iZobs[iZobs]])) * "mm"

							PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], Ψ_Plot[1:N_∑T_Plot, calibr.iZobs[iZobs]],style=Style_Hypix[iZobs], legendentry=Label_Sim, mark="none")
						end # loop

						Plot_ΨhypixTop = map(1:length(CellTop)) do iTop
							Label_HypixTop = "HyPix=" * string(discret.Znode[CellTop[iTop]]) * "mm"

							PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], Ψ_Plot[1:N_∑T_Plot,CellTop[iTop]] ,style=Style_HypixTop[iTop], legendentry=Label_HypixTop, mark="none")
						end # loop

						Ploting_Ψ = [Plot_Ψhypix;]

						push!(Plot_TimeSeries, PGFPlots.Axis(Ploting_Ψ, style="width=20cm, height=8cm", xlabel=L"$\ Time \ [Day]$", ylabel=L"$\Psi \ [cm]$", xmax=∑T_Plot[N_∑T_Plot], ymode="log", legendStyle ="{at={(0.0,-0.3)}, anchor=south west, legend columns=5}"))

				end # if option.hypix.Plot_Ψ


				# PLOTTING GRAPH FLUX
				if option.hypix.Plot_Flux
					Plot_Qhypix = map(1:calibr.Ndepth) do iZobs
						Label_Sim = "HyPix=" * string( (discret.Znode[calibr.iZobs[iZobs]])) * "mm"

						PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], ΔFlux_Plot[1:N_∑T_Plot, calibr.iZobs[iZobs]],style=Style_Hypix[iZobs], legendentry=Label_Sim, mark="none")
					end # loop

					Plot_QhypixTop = map(1:length(CellTop)) do iTop
						Label_HypixTop = "HyPix=" * string(discret.Znode[CellTop[iTop]]) * "mm"

						PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], ΔFlux_Plot[1:N_∑T_Plot,CellTop[iTop]] ,style=Style_HypixTop[iTop], legendentry=Label_HypixTop, mark="none")
					end # loop

					Plot_Recharge = PGFPlots.Plots.Linear(∑T_Plot[1:N_∑T_Plot], -ΔFlux_Plot[1:N_∑T_Plot, N_iZ+1], style="ybar,fill=red, very thick, solid", mark="none", legendentry=L"$\Delta recharge \ [mm \ day^{-1}]$")					

					Ploting_Q = [Plot_Qhypix; Plot_Recharge]

					push!(Plot_TimeSeries, PGFPlots.Axis(Ploting_Q, style="width=20cm, height=8cm", xlabel=L"$\ Time \ [Day]$", ylabel=L"$\Delta Q \ [mm \ day^{-1}]$", xmax=∑T_Plot[N_∑T_Plot], legendStyle ="{at={(0.0,-0.3)}, anchor=south west, legend columns=6}"))

				end # if option.hypix.Plot_Flux

			Path = path.Hypix_calibr * "2_" * string(iSim) * ".svg"	
			PGFPlots.save(Path, Plot_TimeSeries) 

		end  # function: TIME_SERIES


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INTERCEPTION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RAINFALL_INTERCEPTION(clim, ∑T_Climate, iSim)

			Plot_TimeSeries = PGFPlots.GroupPlot(1, 1, groupStyle = "horizontal sep = 3cm, vertical sep = 3cm")

			local ∑T_Int = ceil.(Int, ∑T_Climate[1:clim.N_Climate] .* cst.Second_2_Day)

			Plot_PrThrough = PGFPlots.Plots.Linear(∑T_Int[1:clim.N_Climate], clim.Pr_Through[1:clim.N_Climate], style="const plot, fill=teal, teal", mark="none", legendentry=L"$ Pr_{through} \ [mm \ day^{-1}] $")

			Plot_Pr = PGFPlots.Plots.Linear(∑T_Int[1:clim.N_Climate], clim.Pr[1:clim.N_Climate], style= "const plot, blue", mark="none", legendentry= L"$ Pr \ [mm \ day^{-1}] $")

			Plot_Pet = PGFPlots.Plots.Linear(∑T_Int[1:clim.N_Climate], clim.Pet[1:clim.N_Climate], style= "orange", mark="none", legendentry= L"$ Pet_{int} \ [mm \ day^{-1}] $")

			Plot_Climate = [Plot_Pr, Plot_Pet, Plot_PrThrough]

			push!(Plot_TimeSeries, PGFPlots.Axis(Plot_Climate, style="width=25cm, height=10cm", xlabel=L"$\ Time \ [Day]$", ylabel=L"$Pr \ [cm \ day^{-1}] $",  xmin=0.0, xmax=∑T_Int[clim.N_Climate], ymin=0.0, legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=3}"))

			Path = path.Plot_RainfallInterception * "_" * string(iSim) * ".svg"
			PGFPlots.save(Path, Plot_TimeSeries) 
			
		end  # function: INTERCEPTION


	# ========================================
	# PLOTTING HYDRAULIC RELATIONSHIP FOR EVERY HORIZON
	# ======================================== 
	function θΨK(hydroHorizon, N_iHorizon, iSim)

		# Deriving the Min and Max Ψ from principals of soil physics
		Ψ_Min_Horizon = Vector{Float64}(undef, N_iHorizon)
		Ψ_Max_Horizon = Vector{Float64}(undef, N_iHorizon)
		for iZ=1:N_iHorizon
			Ψ_Max_Horizon[iZ], Ψ_Min_Horizon[iZ] = ΨminΨmax.ΨMINΨMAX(hydroHorizon.θs[iZ], hydroHorizon.θsMacMat[iZ], hydroHorizon.σ[iZ], hydroHorizon.σMac[iZ], hydroHorizon.Ψm[iZ], hydroHorizon.ΨmMac[iZ])
		end  # for iZ=1:N_iHorizon
		
		# PREPARING THE DATA
			N_Se = 1000
			local Ψplot = exp.(range(log(minimum(Ψ_Min_Horizon[1:N_iHorizon])), stop = log(maximum(Ψ_Max_Horizon[1:N_iHorizon])), length=N_Se)) 

			local θplot    = Vector{Float64}(undef, N_Se)
			local Kplot    = Vector{Float64}(undef, N_Se)
			local ∂θ∂Ψplot = Vector{Float64}(undef, N_Se)
			local ∂K∂Ψplot = Vector{Float64}(undef, N_Se)

			Plot_θΨK = PGFPlots.GroupPlot(4, 100, groupStyle = "horizontal sep = 3.5cm, vertical sep = 3.5cm")

		# FOR EVERY HORIZON
		for iZ = 1:N_iHorizon
			
			for iΨ = 1:N_Se
				if Ψ_Max_Horizon[iZ] ≥ Ψplot[iΨ] ≥ Ψ_Min_Horizon[iZ]
					θplot[iΨ]    = wrc.Ψ_2_θDual(Ψplot[iΨ], iZ, hydroHorizon)
					
					Kplot[iΨ]    = kunsat.Ψ_2_KUNSAT(Ψplot[iΨ], iZ, hydroHorizon)
					
					∂θ∂Ψplot[iΨ] = wrc.∂θ∂Ψ(Ψplot[iΨ], iZ, hydroHorizon)

					∂K∂Ψplot[iΨ] = kunsat.∂K∂Ψ(Ψplot[iΨ], iZ, hydroHorizon)
				else
					θplot[iΨ]    = NaN
					
					Kplot[iΨ]    = NaN
					
					∂θ∂Ψplot[iΨ] = NaN

					∂K∂Ψplot[iΨ] = NaN
				end
			end # for iΨ

			Θs_Max = maximum(hydroHorizon.θs[1:N_iHorizon]) + 0.05
			Ks_Min = 10.0 ^ -7 * cst.MmS_2_CmH
			Ks_Max = maximum(hydroHorizon.Ks[1:N_iHorizon]) * cst.MmS_2_CmH * 1.1

			Title =" $(path.SiteName_Hypix)  Horizon = $(iZ)"
		
		# Plot 1: θΨ
			Plot_θΨ = PGFPlots.Plots.Linear(log.(Ψplot) , θplot, style=" smooth, blue, very thick", mark="none", legendentry=L"$ \theta ( \Psi ) $")

			Plot_hydro = [Plot_θΨ]

			push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xlabel=L"$ Ln \ \Psi [mm]$", ylabel=L"$ \theta \ [mm{^3} \ mm^{-3}]$", xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), ymin=0.0, ymax=Θs_Max, legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

		# Plot 2: Kplot(Ψplot)
			Plot_Kθ = PGFPlots.Plots.Linear(log.(Ψplot), Kplot .* cst.MmS_2_CmH, style=" smooth, red, very thick", mark="none", legendentry=L"$ K_{unsat} \ ( \Psi ) $")

			Plot_hydro = [Plot_Kθ]
			
			push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title,  xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), ymax=Ks_Max, ymode="log", xlabel=L"$Ln \  \Psi [mm]$", ylabel=L"$ K_{unsat} \ [cm \ h^{-1}]$", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

		# Plot 3: ∂θ∂Ψplot
			Plot_∂θ∂Ψ = PGFPlots.Plots.Linear(log.(Ψplot), ∂θ∂Ψplot , style=" smooth, green, very thick", mark="none", legendentry=L"$ \partial \theta \partial \Psi $")

			Plot_hydro = [Plot_∂θ∂Ψ]

			push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), xlabel=L"$Ln \  \Psi [mm] $", ylabel=L"$ \partial \theta \partial \Psi $", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

		# Plot 4: ∂K∂Ψplot
			Plot_∂K∂Ψ = PGFPlots.Plots.Linear(log.(Ψplot), ∂K∂Ψplot, style=" smooth, teal, very thick", mark="none", legendentry=L"$ \partial K \partial \Psi $")

			Plot_hydro = [Plot_∂K∂Ψ]

			push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), xlabel=L"$Ln \  \Psi \ [mm]$", ylabel=L"$ \partial K \partial \Psi $", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

		end #iZ ............................................

		Path = path.Plot_Hypix_θΨK * "_" * string(iSim) * ".svg"
		PGFPlots.save(Path, Plot_θΨK) 
	end # function θΨK


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ROOTDENSITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function VEG_FUNCTIONS(discret, iSim, N_iRoot, veg, Z, ΔRootDensity)

			Plot_All = PGFPlots.GroupPlot(2, 1, groupStyle = "horizontal sep = 3cm, vertical sep = 3cm")

			# PLOT VEG_FUNCTIONS
				ΔRootDensity_Norm = Vector{Float64}(undef, N_iRoot)
				# Taking accoung the tickness of the discretisation
				# for iZ=1:N_iRoot
				# 		ΔRootDensity_Norm[iZ] = Z[N_iRoot] * ΔRootDensity[iZ] / discret.ΔZ[iZ]
				# end

				# Plotting
					Plot_RootDensity = PGFPlots.Plots.Linear(ΔRootDensity[1:N_iRoot], discret.Znode[1:N_iRoot], style=" smooth, cyan, very thick", mark="none")

					Plot = [Plot_RootDensity]

					push!(Plot_All, PGFPlots.Axis(Plot, style="width=12cm, height=8cm", xlabel=L"$ \Delta Rdf \ [\%] $", ylabel=L"$Z \ [mm]$", title="(a)"))
			
			# PLOT StressReduction
				# Data	
				N_Se = 6
				Ψstress = Array{Float64}(undef, 2, N_Se) 
				Ψstress[1,1] = veg.Ψfeddes1 / 10.0
				Ψstress[1,2] = veg.Ψfeddes1
				Ψstress[1,3] = veg.Ψfeddes2
				Ψstress[1,4] = veg.Ψfeddes3
				Ψstress[1,5] = veg.Ψfeddes4
				Ψstress[1,6] = veg.Ψfeddes4 * 2.0

				Wsf = Vector{Float64}(undef, N_Se)
				for iΨ ∈ 1:N_Se
					Wsf[iΨ] = rootwateruptake.stressReduction.WATER_STRESS_FUNCTION(2, iΨ, veg, Ψstress)
				end

				Plot_Wsf = PGFPlots.Plots.Linear(Ψstress[1,1:N_Se] .* cst.Mm_2_kPa, Wsf[1:N_Se], style="violet, very thick", mark="none")

				Plot = [Plot_Wsf]

				push!(Plot_All, PGFPlots.Axis(Plot, style="width=12cm, height=8cm", xlabel=L"$ \Psi \ [kPa]$", xmode="log", ylabel=L"$ F_{waterStress} \ [-]$", title="(b)"))

			Path = path.Vegetation * "_" * string(iSim) * ".svg"
			PGFPlots.save(Path, Plot_All)	
		end  # function ROOTDENSITY




			# =============================================================
			#		module: plots
			# =============================================================
			module plots
			import ...sorptivity, ..wrc, ..path, ..cst, ...option, ...param
			export RAINFALL_INTERCEPTION, PLOT_SORPTIVITY

				using Plots.PlotMeasures, LaTeXStrings
				using Plots;pgfplotsx()

				using Dates
				
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				#		FUNCTION : PLOT_SORPTIVITY
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function PLOT_SORPTIVITY(iSim, hydro)
					println("  ==  START: PLOT_SORPTIVITY_SeIni  ==")

					# Setting the range of values for Se
                  Se_Ini         = collect(0.0:0.001:1.0)
                  N_SeIni        = length(Se_Ini)
                  Sorptivity_Mod = Array{Float64}(undef, (N_SeIni))
                  θ_Ini          = Array{Float64}(undef, (N_SeIni))

					for iSeIni=1:N_SeIni
						θ_Ini[iSeIni] = wrc.Se_2_θ(Se_Ini[iSeIni], 1, hydro)

						Sorptivity_Mod[iSeIni] = sorptivity.SORPTIVITY(θ_Ini[iSeIni], 1, hydro) 
					end
					
					# PLOTTING ====================	
						Plot1=Plots.plot(layout=1)

						Title =" $(path.SiteName_Hypix)"

						Plots.plot!(Plot1, Se_Ini[1:N_SeIni] , Sorptivity_Mod[1:N_SeIni], framestyle = [:box :semi :origin :zerolines :grid :true], xlabel=L"Initial \ Se \ [-]", ylabel=L"Sorptivity \  [ \ mm \ \sqrt s \ ]", label="", grid=false) 
					
						Path =path.Plot_Sorptivity  * "_" * string(iSim) * ".svg"

						Plots.savefig(Plot1, Path)

						println("			 ~ ", Path, "~")

				end  # function: PLOT_SORPTIVITY


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# 		FUNCTION : INTERCEPTION
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function RAINFALL_INTERCEPTION(clim, i∑T_CalibrStart_Day, iSim)

					# TICKS
						DateTick=range(clim.Date[i∑T_CalibrStart_Day],step=Day(7),clim.Date[clim.N_Climate])
						
						DateTick2= Dates.format.(DateTick, "d u Y")
					
					# PLOT
						Plot1=Plots.plot(layout=1)

						Title =" $(path.SiteName_Hypix)" 
						
						Plots.plot!(Plot1, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr[i∑T_CalibrStart_Day:clim.N_Climate], color=:blue, colorbar=false,  line = :solid, label= L"$\Delta Pr  $")
						
						Plots.plot!(Plot1, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr_Through[i∑T_CalibrStart_Day:clim.N_Climate], color=:cyan, colorbar=false, label=L"$\Delta Pr_{through}$")

						Plots.plot!(Plot1, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], 10.0*clim.Pet[i∑T_CalibrStart_Day:clim.N_Climate], color=:green, colorbar=false, label= L"$10x\Delta Pet_{int}$")
						
						Plots.plot!(Plot1, grid=false, framestyle=:origin, size=(1000, 600), legend=:topright, xrotation=rad2deg(pi/3), xticks=(DateTick, DateTick2), title=Title, xlabel=L"$Day$", ylabel=L"$Daily \ \Delta Pr  \ \slash \ \Delta Pr_{through} \ \slash \ \Delta Pet_{int} \ [mm] $")

												
					Path = path.Plot_RainfallInterception * "_" * string(iSim) * ".svg"
					Plots.savefig(Plot1, Path)
					println("			 ~ ", Path, "~")
				end  # function: INTERCEPTION

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : TIMESERIES
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function TIMESERIES(∑T_Date_Plot, ∑T_Plot, calibr, discret, Flag_Plot_Pond, iSim, N_∑T_Plot, N_iZ, ΔEvaporation_Plot, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, θ_Plot, θobs_Plot, clim, i∑T_CalibrStart_Day)

				# PATH
					Path = path.Hypix_calibr * "_" * string(iSim) * ".svg"
					rm(Path, force=true, recursive=true)	

				# TICKS
					# Date_Start_Calibr = calibr.Date[1]
					Date_Start_Calibr = DateTime(param.hypix.calibr.Year_Start, param.hypix.calibr.Month_Start, param.hypix.calibr.Day_Start, param.hypix.calibr.Hour_Start, param.hypix.calibr.Minute_Start, param.hypix.calibr.Second_Start) # since we need to compute the culmulativeof the 1rst day
					
					# Date_End_Calibr = calibr.Date[end]
					Date_End_Calibr = DateTime(param.hypix.Year_End, param.hypix.Month_End, param.hypix.Day_End, param.hypix.Hour_End, param.hypix.Minute_End, param.hypix.Second_End)
					
					DateTick=range(Date_Start_Calibr,step=Day(14),Date_End_Calibr)
					
					DateTick2= Dates.format.(DateTick, "d u Y")
				
				# PLOTTING
					Plot = Plots.plot(layout=(3, 1), size=(2500,2200), bottom_margin=0.01mm)
					
					default(titlefont=(20,"times"), legendfontsize=24, guidefont=18, tickfont=18, grid=true)


				# Plot Climate	
				iSubplot = 0
				if option.hypix.Plot_Climate
					iSubplot += 1		

					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], -ΔFlux_Plot[1:N_∑T_Plot, N_iZ+1], label=L"$\Delta Q$", line=(:solid, 1), linecolour=:red, fillcolor=:darkred, fill=(0,:darkred))
					
					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], ΔPond_Plot[1:N_∑T_Plot], label=L"$\Delta H_{Pond}$", linecolour=:grey, fill = (0, :grey))
					
					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr[i∑T_CalibrStart_Day:clim.N_Climate], color=:blue, colorbar=false,  line =(:sticks, :solid, 5), label= L"$\Delta Pr  $")

					Plot_Climate = Plots.plot!(Plot, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr_Through[i∑T_CalibrStart_Day:clim.N_Climate], color=:cyan, line =(:sticks, :solid, 4), colorbar=false, label=L"$\Delta Pr_{through}$")
	
					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ylabel=L"$Daily \ Simulation \ [mm]$", title=path.SiteName_Hypix, xtickfont = (0.01, :white), xrotation=rad2deg(pi/2))
				end # if: option.hypix.Plot_Climate

				# PLOT EVAPOYTRANSPIRATION
					iSubplot += 1	

					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[2:N_∑T_Plot], ΔPet_Plot[2:N_∑T_Plot], linecolour=:darkgreen, label=L"$\Delta Pet$", line=(2.5,:solid))

					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[2:N_∑T_Plot], ΔSink_Plot[2:N_∑T_Plot], linecolour=:red, line=(2.0,:solid), label=L"$\Delta Sink$")

					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[2:N_∑T_Plot], (ΔSink_Plot[2:N_∑T_Plot].-ΔEvaporation_Plot[2:N_∑T_Plot]), label=L"$\Delta Rwu$", linecolour=:blue, line=(2.0,:solid))
					
					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[2:N_∑T_Plot], ΔEvaporation_Plot[2:N_∑T_Plot], label=L"$\Delta Evap$", linecolour=:purple4, line=(2.0,:solid))

					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ylabel=L"$Daily \ Simulation \ [mm]$", xtickfont = (0.01, :white), xrotation=rad2deg(pi/2))

				# PLOT Θ
				if option.hypix.Plot_θ
					iSubplot += 1

					Style_Hypix = [:red, :darkviolet, :orange, :teal, :blue]

					# Observation θplot obs
					for iZobs = 1:calibr.Ndepth
						# lABEL
							Label_Obs = "Obs=" * string(Int(floor(calibr.Z[iZobs]))) * "mm"

							Label_Sim = "Sim=" * string( Int(floor((discret.Znode[calibr.iZobs[iZobs]])))) * "mm"

						# Plotting
							# Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], θobs_Plot[1:N_∑T_Plot, iZobs].+param.hypix.calibr.θobs_Uncert, line=(0.5,:solid), linecolour=Style_Hypix[iZobs], label=false)
		
							# Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], max.(θobs_Plot[1:N_∑T_Plot, iZobs].-param.hypix.calibr.θobs_Uncert, 0.0), line=(0.5,:solid), linecolour=Style_Hypix[iZobs], label=false)

							Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], θobs_Plot[1:N_∑T_Plot, iZobs], line=(2.5,:solid), linecolour=Style_Hypix[iZobs], label=Label_Obs)

							Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], θ_Plot[1:N_∑T_Plot, calibr.iZobs[iZobs]], label=Label_Sim, line=(2.5,:dashdot), linecolour=Style_Hypix[iZobs])
					end # loop

					Plot_θ = Plots.plot!(subplot=iSubplot, ylabel=L"$\theta \ [mm^3 \ mm^{-3}]$")

					Plot = Plots.plot(Plot, Plot_θ, Plot_Climate, xmin=∑T_Date_Plot[1], xmax=∑T_Date_Plot[N_∑T_Plot], ymin=0.0, xtick=(DateTick,DateTick2), xrotation=rad2deg(pi/4), framestyle=:box, grid=true)

				end # if: option.hypix.Plot_θ
				
				Plots.savefig(Plot, Path)
				println("			 ~ ", path.Hypix_calibr, "~")
			
				return nothing
				end  # function: TIMESERIES

			end  # module: plots
			# ............................................................

end  # module plotHypix
# ............................................................