# =============================================================
#		MODULE: hypixStart
# =============================================================

	Path_SoilWater = dirname(dirname(@__DIR__)) * "/src/" # moving down the path twice
	
	include(Path_SoilWater * "Cst.jl")
	include(Path_SoilWater * "Hydro/Wrc.jl")
	include(Path_SoilWater * "Hydro/Kunsat.jl")
	include(Path_SoilWater * "Tool.jl")
	include(Path_SoilWater * "Table.jl")
	include(Path_SoilWater * "Reading.jl")
	include(Path_SoilWater * "Hydro/HydroStruct.jl")
	include(Path_SoilWater * "Hydro/HydroRelation.jl")
	include(Path_SoilWater * "Stats.jl")
	include("HorizonLayer.jl")

	include("ReadWriteHypix/ReadLinkingFile.jl")
	include("ReadWriteHypix/OptionsHypix.jl")
	include("ReadWriteHypix/ParamsHypix.jl")
	include("ReadWriteHypix/PathsHypix.jl")
	include("Opt/ThetaObs.jl")
	include("Climate.jl")
	include("Discretisation.jl")
	include("Memory.jl")
	include("Veg/VegStruct.jl")
	
	include("ReadWriteHypix/ReadHypix.jl")
	include("WaterBalance.jl")
	include("Veg/Interception.jl")
	include("Veg/RootWaterUptake.jl")
	include(Path_SoilWater * "Hydro/ΨminΨmax.jl")
	include("TimeStep.jl")
	include("Evaporation.jl")
	include("Interpolate.jl")
	include("Opt/OfHypix.jl")
	include("Veg/Pet.jl")
	include(Path_SoilWater * "Sorptivity/Sorptivity.jl")
	include("Ponding.jl")
	include("Flux.jl")
	include("Residual.jl")
	include("Richard.jl")
	include("HypixModel.jl")
	include("Opt/HypixOpt.jl")
	include("Other/θaver.jl")
	include("ΔΔtchange.jl")
	# include("θini.jl")

	# include("CheckError.jl")
	# include("Other/PlotHypix.jl")
	# include("Other/PlotOther.jl")

module hypixStart
	import ..cst, ..horizonLayer, ..hydroStruct, ..hypixModel, ..hypixOpt, ..readHypix, ..readLinkingFile, ..stats, ..table, ..waterBalance, ..θaver, ..Δtchange
	import Statistics: mean
	import Dates: now, value

	export HYPIX_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIX_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYPIX_START(ProjectHypix)

			# GETTING PATHS ===

			Time_Start = now()
				Path_Hypix = dirname(dirname(@__DIR__)) # moving down the path twice

				dateHypix, Id, N_Scenario, pathInputHypix, SiteName = readLinkingFile.LINKING_FILE(Path_Hypix, ProjectHypix)

			# READING VALUES FOR EVERY SCENARIOS ===
				paramHypix = []; optionHypix=[]; pathOutputHypix=[]

				for iScenario = 1:N_Scenario #----------------------

					∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑∑ΔSink, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, ∑ΔQ_Bot, CccBest, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Efficiency, Flag_θΨini, Flag_θΨini, Global_WaterBalance, Global_WaterBalance_NormPr, Hpond, hydro, hydro_best, hydroHorizon, hydroHorizon_best, iNonConverge_iOpt, Laiᵀ, Laiᵀ_η, Layer, N_∑T_Climate, N_Layer, NiZ, NseBest, obsTheta, optionHypix, paramHypix, pathInputHypix, pathOutputHypix, Q, Residual, SwcRoots, Temp, veg, veg_best, WilmotBest, WofBest, Z, Zlayer, ΔEvaporation, ΔLnΨmax, ΔPet, ΔPr, ΔRunTimeHypix, ΔSink, ΔT, ΔT_Average, θ, θini_or_Ψini, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest = readHypix.READ_START(dateHypix, Id, iScenario, N_Scenario, Path_Hypix, pathInputHypix, ProjectHypix, SiteName)

				# OPTIMISATION
					if !(optionHypix.opt.Optimisation)
                  paramHypix.iOptMultiStep_Start = 1
                  paramHypix.iOptMultiStep_End   = 1
					end

		for iMultistep = paramHypix.iOptMultiStep_Start:paramHypix.iOptMultiStep_End
			
			println("		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			println("		=== === === START: Looping with time $iMultistep steps \n")

			# COUNT SIMULATIONS
				iOpt_Count = iMultistep - paramHypix.iOptMultiStep_Start + 1

			# OBTAINING HYDRAULIC AND VEGETATION PARAMETERS (depending of we have multistep optimisation)
			if optionHypix.opt.Optimisation
				hydro, hydroHorizon, optim, veg = readHypix.HYPIX_PARAM(Layer, hydro, hydroHorizon, iMultistep, NiZ, optionHypix, paramHypix, pathInputHypix.MultistepOpt[iScenario], veg)

				println("  =============== HyPix Param Read ============== ")

			else
				# options of optim		
					Flag_Opt = false
					NparamOpt = 0

				optim = ( NparamOpt=NparamOpt, Flag_Opt=Flag_Opt)		
			end # optionHypix.Optimisation


			if optim.Flag_Opt
				hydro, hydro_best, hydroHorizon, hydroHorizon_best, veg, veg_best, WofBest = hypixOpt.HYPIXOPTIMISATION_START(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, hydro, hydro_best, hydroHorizon, hydroHorizon_best, iOpt_Count, Laiᵀ, Laiᵀ_η, Layer, N_∑T_Climate, N_Layer, NiZ, obsTheta, optim, optionHypix, paramHypix, Q, Residual, veg, veg_best, WofBest, Z, ΔEvaporation, Hpond, ΔPet, ΔPr, ΔSink, ΔT, ΔLnΨmax, θ, θini_or_Ψini, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest)
			end
			
			# if Flag_Opt then it will rerun with the optimal parameters
			println("\n =============== Start running HyPix =========================== \n")
			∑Pet, ∑Pr, ∑T, ∑T_Climate, clim, discret, iNonConverge, IterCount, N_iRoot, Nit, NiZ, Q, veg, ΔEvaporation, Hpond, ΔRootDensity, ΔT, θ, Ψ = hypixModel.HYPIX_MODEL(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, hydro, Laiᵀ, Laiᵀ_η, N_∑T_Climate, NiZ, optionHypix, paramHypix, Q, Residual, veg, Z, ΔEvaporation, Hpond, ΔPet, ΔPr, ΔSink, ΔT, ΔLnΨmax, θ, θini_or_Ψini, Ψ, Ψ_Max, Ψ_Min, Ψbest)

			# WATER BALANCE
				# Computed after the warmup period
				∑∑WaterBalance, ∑WaterBalance_η, ∑ΔSink, i∑T_CalibrStart, ΔStorage = waterBalance.WATERBALANCE(∑T, obsTheta, discret, hydro, N_iRoot, Nit, NiZ, Q, ΔSink, ΔT, θ, Ψ)

			# SUMMARY HOW GOOD THE SIMULATION
				# Computed climate day after the warmup period
					i∑T_CalibrStart_Day = 1::Int64 
					while ∑T_Climate[i∑T_CalibrStart_Day] < obsTheta.∑T[1]
						i∑T_CalibrStart_Day += 1
					end
					i∑T_CalibrStart_Day += 1

				# Climate
					∑Pr_Clim = sum(clim.Pr[i∑T_CalibrStart_Day:clim.N_Climate]) 
					∑Pet_Net = sum(clim.Pet[i∑T_CalibrStart_Day:clim.N_Climate])

				# Soil water content for the rootzone at the end of simulation
				@fastmath @inbounds @simd for iZ=1:N_iRoot
					SwcRoots[iOpt_Count] += θ[Nit,iZ] * discret.ΔZ[iZ]
				end

				# Timing 
					Time_End = now()
					ΔRunTimeHypix[iOpt_Count] = value(Time_End - Time_Start) / 1000

				# Convergence rate
					iNonConverge_iOpt[iOpt_Count]          = iNonConverge

					Efficiency[iOpt_Count]                 = ceil(Int, cst.Day_2_Second * IterCount / ∑T[Nit] )
					
					ΔT_Average[iOpt_Count]                 = ceil(Int, mean(ΔT[i∑T_CalibrStart:Nit]))
					
					Global_WaterBalance[iOpt_Count]        = ∑∑WaterBalance
					
					Global_WaterBalance_NormPr[iOpt_Count] = 100.0 * ∑WaterBalance_η[Nit]
					
					∑∑ΔSink[iOpt_Count]                    = ∑ΔSink[Nit]
				
				# Ground water recharge
					∑ΔQ_Bot[iOpt_Count] = 0.0
					for iT=i∑T_CalibrStart:Nit
						∑ΔQ_Bot[iOpt_Count] = ∑ΔQ_Bot[iOpt_Count] + ΔT[iT] * Q[iT, NiZ+1]
					end
           
				println("		=== ===START: summary  $iMultistep steps ...")

					println("			∑Pr 			= ", ceil(Int, ∑Pr_Clim), "  [mm]")
					println("			∑Pr_Soil 		= ", ceil(Int, ∑Pr[Nit] - ∑Pr[i∑T_CalibrStart]),  "  [mm]")
					println("			∑Pr_Intercepted/∑Pr 	= ", ceil(Int, 100. * (∑Pr_Clim - (∑Pr[Nit]-∑Pr[i∑T_CalibrStart])) / (∑Pr_Clim + eps(10.0))),  "  [%]")
					println("			∑Pet_Net 		= ", ceil(Int, ∑Pet_Net), "  [mm]")
					println("			∑Pet 			= ", ceil(Int, ∑Pet[Nit]- ∑Pet[i∑T_CalibrStart]), "  [mm]")
					println("			∑ΔSink/∑Pet_Net 	= ", ceil(Int, 100.0 * ∑ΔSink[Nit] /(∑Pet_Net + eps(10.0))), "  [%] \n")
					
					println("			∑SoilWaterContentRootEnd = ", round(SwcRoots[iOpt_Count], digits=3), "  [mm]")
					println("			∑ΔSink 			= ", -ceil(Int, ∑∑ΔSink[iOpt_Count]), "  [mm]")
					println("			∑Infilt_Bot 		= ", -round(∑ΔQ_Bot[iOpt_Count],  digits=5), "  [mm]")
					println("			Hpond at end 		= ", ceil(Int, Hpond[Nit]), "  [mm] \n")

					println("			iNonConverge 			= ", iNonConverge_iOpt[iOpt_Count], "  [count]")
					println("			Global_WaterBalance_NormPr 	= ", round(Global_WaterBalance_NormPr[iOpt_Count], digits=8), "  [%]")
					println("			Global_WaterBalance 		= ", 	round(Global_WaterBalance[iOpt_Count], digits=8), "  [mm]")
					println("			Average ΔT 			= ",  ΔT_Average[iOpt_Count] , "  [seconds]")
					println("			ΔTmin 				= ",   round(minimum(ΔT[i∑T_CalibrStart:Nit]), digits=0) , "  [seconds]")
					println("			ΔTmax 				= ",  round(maximum(ΔT[i∑T_CalibrStart:Nit]), digits=0) , "  [seconds]")
					println("			ΔT_HyPix 			= ", ceil(Int, ΔRunTimeHypix[iOpt_Count]) , "  [seconds]")			
					println("			Efficiency 			= ", Efficiency[iOpt_Count], "  [iTer day-1]")
					println("			Number_of_cells 	     = ", NiZ, "  [-], \n")

					∑T_Reduced, ∑WaterBalanceη_Reduced, Date_Reduced, Nit_Reduced, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔPrGross_Reduced, ΔQ_Reduced, ΔSink_Reduced, ΔT_Reduced, θ_Reduced, θobs_Reduced, Ψ_Reduced = Δtchange.CHANGE_OUTPUT_ΔT(∑Pet[1:Nit], ∑Pr[1:Nit], ∑T[1:Nit], ∑WaterBalance_η[1:Nit], ∑ΔSink[1:Nit], obsTheta, clim, Nit, NiZ, paramHypix, Q[1:Nit,1:NiZ+1], ΔEvaporation[1:Nit], Hpond[1:Nit], ΔT[1:Nit], θ[1:Nit,1:NiZ], Ψ[1:Nit,1:NiZ], ∑T_Climate)

			# Computing average simulated θ to comapre it with average observed θ
			if optionHypix.θavr_RootZone && optionHypix.θobs	
				θsim_Aver = θaver.θAVER(discret; Z=Z, θ_Reduced=θ_Reduced, NiZ=NiZ, Nit_Reduced=Nit_Reduced, Zaver=min(400.0, Z[NiZ]))

			elseif !(optionHypix.θavr_RootZone) && optionHypix.θobs
				# NseBest[iOpt_Count] = stats.NSE(θobs_Reduced[1:Nit_Reduced], θsim_Aver[1:Nit_Reduced])
				# NseBest[iOpt_Count] = ofHypix.θof.RMSE_θ(∑T, obsTheta, Nit, NiZ, θ, θSim)
				θsim_Aver = Float64[]
			
			end

			if  optionHypix.θobs
				for iZobs = 1:obsTheta.Ndepth
					CccBest[iOpt_Count] += stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(θobs_Reduced[1:Nit_Reduced, iZobs], θ_Reduced[1:Nit_Reduced, obsTheta.ithetaObs[iZobs]])

					NseBest[iOpt_Count] += stats.NSE(θobs_Reduced[1:Nit_Reduced, iZobs], θ_Reduced[1:Nit_Reduced, obsTheta.ithetaObs[iZobs]])

					WilmotBest[iOpt_Count] += stats.NSE_WILMOT(θobs_Reduced[1:Nit_Reduced, iZobs], θ_Reduced[1:Nit_Reduced, obsTheta.ithetaObs[iZobs]])
				end # loop
				CccBest[iOpt_Count] = CccBest[iOpt_Count] / obsTheta.Ndepth
				NseBest[iOpt_Count] = NseBest[iOpt_Count] / obsTheta.Ndepth
				WilmotBest[iOpt_Count] = WilmotBest[iOpt_Count] / obsTheta.Ndepth

				println("			CccBest 			= ", round(CccBest[iOpt_Count], digits=5))
				println("			NseBest 			= ", round(NseBest[iOpt_Count], digits=5))
				println("			WilmotBest 	    = ", round(WilmotBest[iOpt_Count], digits=5))
			end	

			println("		=== === END: summary \n")


			if optionHypix.Table
			println("		=== === START: Table === ===")
	
				# Writing values of hydraulic parameters
				table.hyPix.HYDRO(hydroHorizon, iMultistep, N_Layer, pathOutputHypix)

				# Writing values of veg parameters
				table.hyPix.VEG(veg, iMultistep, pathOutputHypix)

				table.hyPix.PERFORMANCE(∑∑ΔSink, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iOpt, iMultistep, iScenario, NseBest, paramHypix,  pathOutputHypix, SwcRoots, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average)	
					

				if optionHypix.Table_Discretization
					table.hyPix.DISCRETISATION_RRE(discret, NiZ, Z[1:NiZ], pathOutputHypix)
				end
				if  optionHypix.Table_TimeSeries
					table.hyPix.θDATA(∑T[1:Nit], ΔT[1:Nit], ∑Pr[1:Nit], ΔPr[1:Nit], Hpond[1:Nit], ΔT[1:Nit].*Q[1:Nit,1], ∑WaterBalance_η[1:Nit], iMultistep, pathOutputHypix)
				end
				if optionHypix.Table_TimeSeriesDaily
					table.hyPix.TIME_SERIES_DAILY(∑T_Reduced[1:Nit_Reduced], ∑WaterBalanceη_Reduced[1:Nit_Reduced], Date_Reduced[1:Nit_Reduced], iMultistep, Nit_Reduced, ΔEvaporation_Reduced[1:Nit_Reduced], ΔQ_Reduced[1:Nit_Reduced, NiZ+1], ΔPet_Reduced[1:Nit_Reduced], ΔPond_Reduced[1:Nit_Reduced], ΔPr_Reduced[1:Nit_Reduced], ΔSink_Reduced[1:Nit_Reduced], pathOutputHypix)
				end
				if optionHypix.Table_θ
					table.hyPix.θ(Date_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced,1:NiZ], discret.Znode[1:NiZ], iMultistep, pathOutputHypix)
				end
				if optionHypix.Table_Ψ
					table.hyPix.Ψ(Date_Reduced[1:Nit_Reduced], Ψ_Reduced[1:Nit_Reduced,1:NiZ], discret.Znode[1:NiZ], iMultistep, pathOutputHypix)
				end
				if optionHypix.Table_Q
					table.hyPix.Q(Date_Reduced[1:Nit_Reduced], ΔQ_Reduced[1:Nit_Reduced,1:NiZ+1], Z[NiZ], discret.Znode[1:NiZ], iMultistep, pathOutputHypix)
				end
				if optionHypix.Tabule_θΨ
					table.hyPix.θΨ(hydroHorizon, iMultistep, N_Layer, optionHypix, paramHypix, pathOutputHypix)
					table.hyPix.KΨ(hydroHorizon, iMultistep, N_Layer, optionHypix, paramHypix, pathOutputHypix)
				end
				if optionHypix.Table_Climate
					table.hyPix.DAILY_CLIMATE(∑T_Climate, clim, iMultistep, pathOutputHypix)
				end
				if optionHypix.θavr_RootZone && optionHypix.θobs
					table.hyPix.θAVERAGE(Date_Reduced[1:Nit_Reduced], iMultistep, θobs_Reduced[1:Nit_Reduced], θsim_Aver[1:Nit_Reduced], pathOutputHypix)
				end
			println("		=== === END: Table === === \n")
			end  # if optionHypix.Table

		
			if optionHypix.Ploting
			println("		=== === START: Plotting === ===")

				# if optionHypix.Plot_Other
				
				# plotOther.plots.WOF_STEPS(path)
				# plotOther.PLOT_θΨ_Δθ(hydro, pathOutputHypix, paramHypix, optionHypix)
				# 	# plotOther.ΨMINΨMAX(hydro, pathOutputHypix)
				
				# 	# plotOther.SE_Ψ_CONSTRAINED(hydro, pathOutputHypix)
				# 	# plotOther.PLOT_σ_2_θr(hydro, pathOutputHypix)
				# 	# plotOther.PLOT_θΨ_Δθ(hydro, pathOutputHypix)
				# end # optionHypix.Plot_Other

				if optionHypix.Plot_Hypix
					plotHypix.makkie.TIMESERIES(∑T_Reduced, clim, Date_Reduced, discret, i∑T_CalibrStart_Day, iMultistep, iScenario, Nit_Reduced, NiZ, obsTheta, optionHypix, paramHypix,  pathOutputHypix, SiteName, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔPrGross_Reduced, ΔQ_Reduced, ΔSink_Reduced, θ_Reduced, θobs_Reduced, θsim_Aver)

					if optionHypix.Plot_θprofile
						plotHypix.makkie.θPROFILE(∑T_Reduced, discret, iScenario, NiZ, obsTheta, optionHypix, paramHypix, pathOutputHypix, SiteName, θ_Reduced)
					end  # if: optionHypix.Plot_
				end

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
			println("		=== === END: Plotting === === \n")
			end # if optionHypix.Plotting
	
			println("	=== === === END  : Looping with time ")
			println("	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n")
		end # for loop: iMultistep

		println("\n ==== END RUNNING Hypix_1D =====")
			end # for iScenario = 1:N_Scenario
				
		end  # function: HYPIX_START
	# ------------------------------------------------------------------



end  # module hydro
# ............................................................

@time hypixStart.HYPIX_START("LYSIMETERS")