# =============================================================
#		module: plotOther
# =============================================================
module plotOther
	

	# =============================================================
	#		module: plots
	# =============================================================
	module plots
		import ...tool
		using Plots, Plots.PlotMeasures, LaTeXStrings
		using Plots; pgfplotsx()
		# using PGFPlots

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : WOF_STEPS
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function WOF_STEPS(path)

				Path_Output = pathHypix.Plot_OfStep * "Multiplestep.svg"

				rm(Path_Output, force=true, recursive=true)	

				Label = ["Waitoa";"Waihou";"Taupo";"Otorohanga";"Hamilton"]

				Plot1=Plots.plot(layout=(2,2), size=(1000,600), bottom_margin=20px, grid=:x)

				# TICKS 
					Ticks =[ "1","2a","2b","3a","3b","4a","4b","5a","5b"]

				# OF STEP
					Path = pathHypix.Input_OfStep * "Of_Step.CSV"

					Of_Waitoa, N_Waitoa = tool.readWrite.READ_HEADER(Path, "Waitoa")
					Of_Waihou, N_Waihou = tool.readWrite.READ_HEADER(Path, "Waihou")
					Of_Taupo, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Taupo")
					Of_Otorohanga, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Otorohanga")
					Of_Hamilton, N_Hamilton = tool.readWrite.READ_HEADER(Path, "Hamilton")
					
					Id=1:N_Waitoa

					Of = [Of_Waitoa Of_Waihou Of_Taupo Of_Otorohanga Of_Hamilton]

					for i=1:5
						Plots.plot!(Plot1, subplot=1, Id, Of[1:N_Waitoa, i], palette=:darkrainbow,  marker=(:circle, 4, 1.0), line=(2.0,:solid))
					end
					Plots.plot!(Plot1, subplot=1,  xlabel="", ylabel=L"WOF _{\theta} \ [mm^{3} \ mm^{-3}]", xticks=(1:1:9, Ticks), xtickfont=(12, :white), legend=false, title="(a) Weighted Objective Function", titlelocation = :left)

				# GROUNDWATER STEP
					Path = pathHypix.Input_OfStep * "Groundwater_Step.csv"

					Groundwater_Waitoa, N_Waitoa = tool.readWrite.READ_HEADER(Path, "Waitoa")
					Groundwater_Waihou, N_Waihou = tool.readWrite.READ_HEADER(Path, "Waihou")
					Groundwater_Taupo, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Taupo")
					Groundwater_Otorohanga, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Otorohanga")
					Groundwater_Hamilton, N_Hamilton = tool.readWrite.READ_HEADER(Path, "Hamilton")

					Groundwater = [Groundwater_Waitoa Groundwater_Waihou Groundwater_Taupo Groundwater_Otorohanga Groundwater_Hamilton]

					for i=1:5
						Plots.plot!(Plot1, subplot=2, Id, Groundwater[1:N_Waitoa, i], palette=:darkrainbow,  marker=(:circle, 4, 1.0), line=(2.0,:solid))
					end
					# Plots.plot!(Plot1, subplot=2, xlabel=L"Multistep \ Optimisation \ Steps", ylabel=L"\zeta _{Q} ", xticks=(1:1:8, Ticks), legend=false, title="Groundwater", titlelocation = :left)

					Plots.plot!(Plot1, subplot=2,  xlabel="", ylabel=L"\zeta _{Q} \ [\%]", xticks=(1:1:9, Ticks), xtickfont=(12, :white), legend=false, title="(b) Drainage", titlelocation = :left)


				# EvapoTranspiration STEP
					Path = pathHypix.Input_OfStep * "Sink_Step.csv"

					Sink_Waitoa, N_Waitoa = tool.readWrite.READ_HEADER(Path, "Waitoa")
					Sink_Waihou, N_Waihou = tool.readWrite.READ_HEADER(Path, "Waihou")
					Sink_Taupo, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Taupo")
					Sink_Otorohanga, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Otorohanga")
					Sink_Hamilton, N_Hamilton = tool.readWrite.READ_HEADER(Path, "Hamilton")

					Sink = [Sink_Waitoa Sink_Waihou Sink_Taupo Sink_Otorohanga Sink_Hamilton]

					for i=1:5
						Plots.plot!(Plot1, subplot=3, Id, Sink[1:N_Waitoa, i],  palette=:darkrainbow,  marker=(:circle, 4, 1.0), line=(2.0,:solid))
					end
					# Plots.plot!(Plot1, subplot=3, xlabel=L"Multistep \ Optimisation \ Steps?", ylabel=L"\zeta _{et} ", xticks=(1:1:8, Ticks), legend=false, title="EvaopoTranspiration", titlelocation = :left)

					Plots.plot!(Plot1, subplot=3, xlabel=L"Multistep \ optimization \ [Layers]", ylabel=L"\zeta _{et}  \ [\%]", xticks=(1:1:9, Ticks), tickfont=(12, :black), legend=false, title="(c) Evapotranspiration", titlelocation = :left)

				# Soil Water Content STEP
					Path = pathHypix.Input_OfStep * "Swc_Step.csv"

					Swc_Waitoa, N_Waitoa = tool.readWrite.READ_HEADER(Path, "Waitoa")
					Swc_Waihou, N_Waihou = tool.readWrite.READ_HEADER(Path, "Waihou")
					Swc_Taupo, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Taupo")
					Swc_Otorohanga, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Otorohanga")
					Swc_Hamilton, N_Hamilton = tool.readWrite.READ_HEADER(Path, "Hamilton")

					Swc = [Swc_Waitoa Swc_Waihou Swc_Taupo Swc_Otorohanga Swc_Hamilton]

					for i=1:5
						Plots.plot!(Plot1, subplot=4, Id, Swc[1:N_Waitoa, i], palette=:darkrainbow, marker=(:circle, 4, 1.0), label=Label[i], line=(2.0,:solid))
					end
					# Plots.plot!(Plot1, subplot=4, xlabel=L"Multistep \ Optimisation \ Steps", ylabel=L"\zeta _{\theta} ", xticks=(1:1:8, Ticks), legend=(-0.15,-0.18), title="Soil Water Content", titlelocation = :left)

					Plots.plot!(Plot1, subplot=4, xlabel=L"Multistep \ optimization \ [Layers]", ylabel=L"\zeta_{swc}  \ [\%]", xticks=(1:1:9, Ticks), tickfont=(12, :black), legend=(0.75,1.0), title="(d) Root zone soil water content", titlelocation = :left)

					Plots.savefig(Plot1, Path_Output)
					println("			 ~ ", Path_Output, "~")
			end  # function: WGroundwater_STEPS


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ROOTDENSITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function VEG_FUNCTIONS(discret, iMultistep, N_iRoot, veg, Z, ΔRootDensity, pathHyPix)

		# 	Plot_All = PGFPlots.GroupPlot(2, 1, groupStyle = "horizontal sep = 3cm, vertical sep = 3cm")

		# 	# PLOT VEG_FUNCTIONS
		# 		ΔRootDensity_Norm = fill(0.0::Float64, N_iRoot)
		# 		# Taking accoung the tickness of the discretisation
		# 		# for iZ=1:N_iRoot
		# 		# 		ΔRootDensity_Norm[iZ] = Z[N_iRoot] * ΔRootDensity[iZ] / discret.ΔZ[iZ]
		# 		# end

		# 		# Plotting
		# 			Plot_RootDensity = PGFPlots.Plots.Linear(ΔRootDensity[1:N_iRoot], discret.Znode[1:N_iRoot], style=" smooth, green, very thick", mark="none")

		# 			Plot = [Plot_RootDensity]

		# 			push!(Plot_All, PGFPlots.Axis(Plot, style="width=12cm, height=8cm", xlabel=L"$ \Delta Rdf \ [\%] $", ylabel=L"$Z \ [mm]$", title="(a)"))
			
		# 	# PLOT StressReduction
		# 		# Data	
		# 		N_Se = 6
		# 		Ψstress = fill(0.0::Float64, 2, N_Se) 
		# 		Ψstress[1,1] = -veg.Ψfeddes1 / 10.0
		# 		Ψstress[1,2] = -veg.Ψfeddes1
		# 		Ψstress[1,3] = -veg.Ψfeddes2
		# 		Ψstress[1,4] = -veg.Ψfeddes3
		# 		Ψstress[1,5] = -veg.Ψfeddes4
		# 		Ψstress[1,6] = -veg.Ψfeddes4 * 2.0

		# 		Wsf = fill(0.0::Float64, N_Se)
		# 		for iΨ ∈ 1:N_Se
		# 			Wsf[iΨ] = rootWaterUptake.stressReduction.WATER_STRESS_FUNCTION(2, iΨ, veg, Ψstress)
		# 		end

		# 		Plot_Wsf = PGFPlots.Plots.Linear(-Ψstress[1,1:N_Se] .* cst.Mm_2_kPa, Wsf[1:N_Se], style="violet, very thick", mark="none")

		# 		Plot = [Plot_Wsf]

		# 		push!(Plot_All, PGFPlots.Axis(Plot, style="width=12cm, height=8cm", xlabel=L"$ \Psi \ [kPa]$", xmode="log", ylabel=L"$ F_{waterStress} \ [-]$", title="(b)"))

		# 	Path = pathHyPix.Vegetation * "_" * string(iMultistep) * ".svg"
		# 	PGFPlots.save(Path, Plot_All)	
		# end  # function ROOTDENSITY
		
	end  # module: plots
	# ............................................................



	

	# =============================================================
	#		module: makkie
	# =============================================================
	module makkie

	import ...hydroStruct, ...hydroRelation, ...wrc, ...timeStep
	using CairoMakie
	export PLOT_θΨ_Δ, DRY_METHOD

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT_θψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_θΨ_Δθ(hydro, pathHyPix, paramHypix, optionHypix)
			N_Se = 100

			optionₘ = optionHypix

			hydroHorizon₂ = hydroStruct.HYDROSTRUCT(optionₘ, 1)
			
			hydroHorizon₂.θs[1] = 0.5
			hydroHorizon₂.θsMacMat[1] = 0.5
			hydroHorizon₂.θr[1] = 0.00
			hydroHorizon₂.σ[1] = 2.0
			hydroHorizon₂.Ψm[1] = hydroRelation.σ_2_Ψm(hydroHorizon₂.σ[1], paramHypix.hydro.kg.Ψσ,  hydro.Ψm_Min[1],  hydro.Ψm_Max[1])
			hydroHorizon₂.ΨmMac[1] = 100.0
			hydroHorizon₂.σMac[1] = 2.0

			# Feasible range 
				Ψ_Max= exp(16.0)
				Ψ_Min= 0.0
			# Plotting the curve
				Ψplot = exp.(range(log(Ψ_Min), stop=log(Ψ_Max), length=N_Se)) 
				θplot = fill(0.0::Float64, N_Se)
				for iΨ = 1:N_Se
					θplot[iΨ] = wrc.Ψ_2_θDual(optionₘ, Ψplot[iΨ], 1, hydroHorizon₂)			
				end # for iΨ

			# INITIALIZING PLOT
				CairoMakie.activate!()
				Makie.inline!(true)

				# AX1
	
				Fig = Figure(resolution = (2000, 800))
				
				Title = L"Time-stepping: $\Delta \theta$"
					Ax1 = Axis(Fig[1,1], title=Title, xlabel= L"ln $\psi$ $[mm]$", ylabel=  L"$\theta$ $[m^{3}  m^{-3}]$",  font="Computer Modern", titlesize=30, fontsize=30, xlabelsize=30, ylabelsize=30 , xgridvisible=false, ygridvisible=false, ytickalign=0.1)

					noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
					Label(Fig[1, 1, TopLeft()], "    (A)", textsize = 24, font = noto_sans_bold, padding = (0, 5, 5, 0), halign = :right)

					# xlims!(Ax1, log(Ψ_Min), log(Ψ_Max))
					xlims!(Ax1, -log(Ψ_Max), -log(Ψ_Min))
					ylims!(Ax1, 0.0, hydroHorizon₂.θs[1]+0.1)
					Ax1.yticks=0:0.1:hydroHorizon₂.θs[1]

					lines!(Ax1, -log.(Ψplot), θplot, linewidth=2, color=:red)
					# Plotting points on the curve
						N_Δθ = 6
						# Ψ_Δθ = fill(0.0::Float64, N_Δθ)
						# Δθ = range(hydroHorizon₂.θr[1], hydroHorizon₂.θs[1] , length=N_Δθ)

						θhalf = (hydroHorizon₂.θs[1] - hydroHorizon₂.θr[1]) / 2.0
						Δθ= 1.0E-2

7
						Δθ = [θhalf-Δθ, θhalf+Δθ, hydroHorizon₂.θs[1], hydroHorizon₂.θs[1]-Δθ, hydroHorizon₂.θr[1], hydroHorizon₂.θr[1]+Δθ]

						for iθ = 1:N_Δθ
							Ψ_Δθ = wrc.θ_2_ΨDual(optionₘ, Δθ[iθ], 1, hydroHorizon₂)

							scatter!(Ax1, [-log(Ψ_Δθ),-log(Ψ_Δθ)], [Δθ[iθ],Δθ[iθ]], markersize=10, color=:blue)
							# scatter!(Ax1, [log(hydroHorizon₂.Ψm[1])], [0.5 * (hydroHorizon₂.θsMacMat[1]+hydroHorizon₂.θr[1])], markersize=20,  color=:red)
							lines!(Ax1, [-log(Ψ_Δθ), -log(Ψ_Δθ)], [0.0, Δθ[iθ]], linestyle=:dash, linewidth=2, color=:grey)
							lines!(Ax1, [-log(Ψ_Δθ), -log(Ψ_Min)], [Δθ[iθ], Δθ[iθ]], linestyle=:dash, linewidth=2, color=:grey)
						end

				# AX2 ============================================
				Title = L"Time-stepping: $\Delta \Psi$"

				Ax2 = Axis(Fig[1,2], title=Title, xlabel= L"ln $\psi$ $[mm]$",  font="Computer Modern", titlesize=25, fontsize=16, xlabelsize=22, ylabelsize=22 , xgridvisible=false, ygridvisible=false,  ytickalign=0.1)

				noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
				Label(Fig[1, 2, TopLeft()], "(B)", textsize = 24, font = noto_sans_bold, padding = (5, 5, 5, 5), halign = :right)

				xlims!(Ax2, -log(Ψ_Max), -log(Ψ_Min))
				ylims!(Ax2, 0.0, hydroHorizon₂.θs[1])
				Ax2.yticks=0:0.1:hydroHorizon₂.θs[1]

				hideydecorations!(Ax2, ticks=false, grid=false)

				lines!(Ax2, -log.(Ψplot), θplot, linewidth=2, color=:red)

				# ΔLnΨmax = fill(0.0: 1)
				# ΔLnΨmax = timeStep.ΔΨMAX!(hydroHorizon₂, 1, optionHypix, paramHypix, ΔLnΨmax)

				# timeStep.ΔθMAX(hydro, 1, 1, optionHypix, ΔLnΨmax, Ψ)
				# Plotting points on the curve
					N_Ψ = 5
					Δθ = fill(0.0::Float64, N_Ψ)

					ΔLogΨ = range(2.0, 16.0 , length=N_Ψ)


					for iΨ = 1:N_Ψ
						Ψ_Δθ = exp(ΔLogΨ[iΨ])

						Δθ[iΨ] = wrc.Ψ_2_θDual(optionₘ, Ψ_Δθ, 1, hydroHorizon₂)

						scatter!(Ax2, [-log(Ψ_Δθ),-log(Ψ_Δθ)], [Δθ[iΨ],Δθ[iΨ]], markersize=10,  color=:blue)
						# scatter!(Ax2, [log(hydroHorizon₂.Ψm[1])], [0.5 * (hydroHorizon₂.θsMacMat[1]+hydroHorizon₂.θr[1])], marker = "x", markersize=15,  color=:red)
						lines!(Ax2, [-log(Ψ_Δθ), -log(Ψ_Δθ)], [0.0, Δθ[iΨ]], linestyle=:dash, linewidth=2, color=:grey)
						lines!(Ax2, [-log(Ψ_Δθ), -log(Ψ_Min)], [Δθ[iΨ], Δθ[iΨ]], linestyle=:dash, linewidth=2, color=:grey)
					end

				trim!(Fig.layout)
				colgap!(Fig.layout, 1)
				rowgap!(Fig.layout, 1)

				save( pathHyPix.Plot_θΨ_Δθ, Fig)
				println("			 ~ ", pathHyPix.Plot_θΨ_Δθ, "~")
				display(Fig)

		return nothing
		end  # function: PLOT_θψ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : name
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DRY_METHOD()

		σobs=	[0.767528364,
					0.84759782,
					1.154454689,
					1.56910817,
					1.959401769,
					1.850609728,
					1.554592632,
					2.084719856,
					2.436491598,
					2.416681323,
					3.24090118,
					3.253171388]

		Ψwet_Obs =	[26,
					23,
					23,
					21,
					17,
					18,
					9,
					8,
					5,
					2,
					0.01,
					0.01]	

		Ψdry_Obs = [2090,
					3550,
					8010,
					18650,
					35270,
					29420,
					15770,
					28230,
					34260,
					24410,
					19490,
					20750]

			N = 100
			σmat = range(0.8, 4, length=N)

			Ψwet = fill(0.0, N)
			Ψdry =  fill(0.0, N)

			for  (i , σ) in enumerate(σmat)
				Ψwet[i] =  max(-2.3116 * σ ^ 2.0 - 2.9372 * σ + 27.83, 0.0)

				Ψdry[i] = exp(1.6216 * log(σ) + 8.7268)
			end

			# INITIALIZING PLOT
				CairoMakie.activate!()
				Makie.inline!(true)


				Fig = Figure(resolution = (1000, 600))

				Ax1 = Axis(Fig[1,1], xlabel= L"$\sigma$ $[-]$", ylabel=  L"$\Psi$ $[mm]$",  font="Computer Modern", titlesize=35, fontsize=35, xlabelsize=35, ylabelsize=35 , xgridvisible=false, ygridvisible=false, yscale = log10, yminorticksvisible = true, yminorgridvisible = true,
				yminorticks = IntervalsBetween(10))

				Ax1.yticks=[0, -1, -100,-1000,-10000,-100000]
				Ax1.yscale = Makie.pseudolog10

				lines!(Ax1, σmat, -Ψwet, linewidth=3, color=:blue, label = L"\Psi _{wet} Model")
				lines!(Ax1, σmat, -Ψdry, linewidth=3, color=:red, label = L"\Psi _{dry} Model")
				scatter!(Ax1, σobs, -Ψwet_Obs, markersize=10, color=:blue, label = L"\Psi _{wet} Zhu")
				scatter!(Ax1, σobs, -Ψdry_Obs, markersize=10, color=:red, label = L"\Psi _{dry} Zhu")

				Leg = Legend(Fig[2,1], Ax1, framevisible=true, orientation=:horizontal, tellheight=true, nbanks=1, framecolor = (:grey,0.5), labelsize=20)

				trim!(Fig.layout)

				display(Fig)
				save( "D:/Main/MODELS/SoilWater_ToolBox/data/OUTPUT/Hypix/RESULTS/Fig3 Zhu_DryMethod.svg", Fig)
		return nothing
		end  # function: name
		# ------------------------------------------------------------------
		
	end  # module: makkie
	# ............................................................

	
	
end  # module: plotOther
# ............................................................