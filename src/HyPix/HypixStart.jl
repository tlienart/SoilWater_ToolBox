# =============================================================
#		MODULE: hypixStart
# =============================================================

include("IncludeHypix.jl")

module hypixStart

	import ..cst, ..horizonLayer, ..hydroStruct, ..hypixModel, ..hypixOpt, ..readHypix, ..readLinkingFile, ..stats, ..tableHypix, ..waterBalance, ..θaver, ..Δtchange
	import Statistics: mean
	import Dates: now, value
	import ..plotHypix, ..tool

	export HYPIX_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIX_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYPIX_START(ProjectHypix)

		println("\n =============== Start running HyPix 1D =========================== \n")

		Time_Start = now()
		# GETTING PATHS ===
		
			Path_Hypix = dirname(dirname(@__DIR__)) # moving down the path twice

			dateHypix, Id, N_Scenario, pathInputHypix, SiteName = readLinkingFile.LINKING_FILE(Path_Hypix, ProjectHypix)

		# READING VALUES FOR EVERY SCENARIOS ===
			for iScenario = 1:N_Scenario #----------------------

				println("=== === === SITENAME=  ", ProjectHypix,"_", SiteName[iScenario], "  === === ===  \n")

				∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑∑ΔSink, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, ∑ΔQ_Bot, CccBest, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Efficiency, Flag_θΨini, Flag_θΨini, Global_WaterBalance, Global_WaterBalance_NormPr, Hpond, hydro, hydro_best, hydroHorizon, hydroHorizon_best, iNonConverge_iOpt, Laiᵀ, Laiᵀ_η, Layer, N_∑T_Climate, N_Layer, NiZ, NseBest, obsTheta, optionHypix, paramHypix, pathInputHypix, pathOutputHypix, Q, Residual, SwcRoots, Temp, veg, veg_best, WilmotBest, WofBest, Z, Zlayer, ΔEvaporation, ΔLnΨmax, ΔPet, ΔPr, ΔRunoff, ΔRunTimeHypix, ΔSink, ΔT, ΔT_Average, θ, θini_or_Ψini, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest = readHypix.READ_START(dateHypix, Id, iScenario, N_Scenario, Path_Hypix, pathInputHypix, ProjectHypix, SiteName)

				# OPTIMISATION
					if !(optionHypix.opt.Optimisation)
						paramHypix.opt.iOptMultiStep_Start = 1
						paramHypix.opt.iOptMultiStep_End   = 1
					end

				for iMultistep = paramHypix.opt.iOptMultiStep_Start:paramHypix.opt.iOptMultiStep_End
					
				if optionHypix.opt.Optimisation
					println("		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
					println("		=== === === START: Multistep= ", iMultistep, "steps, \n")
				end

					# COUNT SIMULATIONS
						iOpt_Count = iMultistep - paramHypix.opt.iOptMultiStep_Start + 1

					# OBTAINING HYDRAULIC AND VEGETATION PARAMETERS (depending of we have multistep optimisation)
					if optionHypix.opt.Optimisation
						hydro, hydroHorizon, optim, veg = readHypix.HYPIX_PARAM_OPT(Layer, hydro, hydroHorizon, iMultistep, NiZ, optionHypix, paramHypix, pathInputHypix.MultistepOpt[iScenario], veg)

						# INJECTING HYDRAULIC & VEGETATION PARAMETERS FROM OUTPUT
						if optionHypix.opt.HydroVegParamReadFromOutput && paramHypix.opt.iOptMultiStep_Start == iMultistep
							@warn( "\n ========= HydroVegParamReadFromOutput = true ========= \n" )

							# Reading hydraulic parameters
								Path = pathOutputHypix.Table_Hydro  * "_" * string(paramHypix.opt.iOptMultiStep_Start-1) * ".csv"
								hydroHorizon, ~ = tool.readWrite.READ_STRUCT_SIMPLE(hydroHorizon, Path)
								@inbounds @simd for iZ=1:N_Layer
									hydroHorizon.So[iZ] = paramHypix.So # 1.0E-8
								end
								hydro = horizonLayer.HYDROHORIZON_2_HYDRO(hydroHorizon, Layer, NiZ, optionHypix)

							# Reading vegetation parameters
								Path = pathOutputHypix.Table_Veg * "_" * string(paramHypix.opt.iOptMultiStep_Start-1) * ".csv"
								veg, ~ = tool.readWrite.READ_STRUCT_SIMPLE(veg, Path)
						end # INJECTING HYDRAULIC & VEGETATION PARAMETERS FROM OUTPUT

					else
						# options of optim		
							Flag_Opt = false
							NparamOpt = 0

						optim = ( NparamOpt=NparamOpt, Flag_Opt=Flag_Opt)		
					end # optionHypix.Optimisation


					if optim.Flag_Opt
						hydro, hydro_best, hydroHorizon, hydroHorizon_best, veg, veg_best, WofBest = hypixOpt.HYPIXOPTIMISATION_START(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, hydro, hydro_best, hydroHorizon, hydroHorizon_best, iOpt_Count, Laiᵀ, Laiᵀ_η, Layer, N_∑T_Climate, N_Layer, NiZ, obsTheta, optim, optionHypix, paramHypix, Q, Residual, veg, veg_best, WofBest, Z, ΔEvaporation, Hpond, ΔPet, ΔPr, ΔRunoff, ΔSink, ΔT, ΔLnΨmax, θ, θini_or_Ψini, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest)
					end
				
					# if Flag_Opt then it will rerun with the optimal parameters
					∑Pet, ∑Pr, ∑T, ∑T_Climate, clim, discret, Hpond, iNonConverge, IterCount, N_iRoot, Nit, NiZ, Q, veg, ΔEvaporation, ΔRootDensity, ΔRunoff, ΔT, θ, Ψ = hypixModel.HYPIX_MODEL(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, hydro, Laiᵀ, Laiᵀ_η, N_∑T_Climate, NiZ, optionHypix, paramHypix, Q, Residual, veg, Z, ΔEvaporation, Hpond, ΔPet, ΔPr, ΔRunoff, ΔSink, ΔT, ΔLnΨmax, θ, θini_or_Ψini, Ψ, Ψ_Max, Ψ_Min, Ψbest)

					# WATER BALANCE
					# Computed after the warmup period
						∑∑WaterBalance, ∑WaterBalance_η, ∑ΔSink, i∑T_CalibrStart, ΔStorage = waterBalance.WATERBALANCE(∑T, obsTheta, discret, hydro, N_iRoot, Nit, NiZ, Q, ΔSink, ΔT, θ, Ψ)

						∑Runoff_Reduced, ∑T_Reduced, ∑WaterBalanceη_Reduced, Date_Reduced, Nit_Reduced, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔPrGross_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔSink_Reduced, ΔT_Reduced, θ_Reduced, θobs_Reduced, Ψ_Reduced = Δtchange.CHANGE_OUTPUT_ΔT(∑Pet[1:Nit], ∑Pr[1:Nit], ∑T[1:Nit], ∑WaterBalance_η[1:Nit], ∑ΔSink[1:Nit], obsTheta, clim, Nit, NiZ, paramHypix, Q[1:Nit,1:NiZ+1], ΔEvaporation[1:Nit], ΔRunoff[1:Nit], Hpond[1:Nit], ΔT[1:Nit], θ[1:Nit,1:NiZ], Ψ[1:Nit,1:NiZ], ∑T_Climate)

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
							∑ΔQ_Bot[iOpt_Count] = 0.0::Float64
							for iT=i∑T_CalibrStart:Nit
								∑ΔQ_Bot[iOpt_Count] = ∑ΔQ_Bot[iOpt_Count] + ΔT[iT] * Q[iT, NiZ+1]
							end
				
						println("		=== ===START: summary  $iMultistep steps ...")

						if optionHypix.TopBoundary⍰ ≠ "Ψ"
							println("			∑Pr 			= ", ceil(Int, ∑Pr_Clim), "  [mm]")
							println("			∑Pr_Soil 		= ", ceil(Int, ∑Pr[Nit] - ∑Pr[i∑T_CalibrStart]),  "  [mm]")
							println("			∑Pr_Intercepted/∑Pr 	= ", ceil(Int, 100. * (∑Pr_Clim - (∑Pr[Nit]-∑Pr[i∑T_CalibrStart])) / (∑Pr_Clim + eps())),  "  [%]")
						end
						if optionHypix.RootWaterUptake
							println("			∑Pet_Net 		= ", ceil(Int, ∑Pet_Net), "  [mm]")
							println("			∑Pet 			= ", ceil(Int, ∑Pet[Nit]- ∑Pet[i∑T_CalibrStart]), "  [mm]")
							println("			∑ΔSink/∑Pet_Net 	= ", ceil(Int, 100.0 * ∑ΔSink[Nit] /(∑Pet_Net + eps(10.0))), "  [%] \n")
						end
						
						println("			∑SoilWaterContentRootEnd = ", round(SwcRoots[iOpt_Count], digits=3), "  [mm]")
						println("			∑ΔSink 			= ", -ceil(Int, ∑∑ΔSink[iOpt_Count]), "  [mm]")
						println("			∑Infilt_Bot 		= ", -round(∑ΔQ_Bot[iOpt_Count],  digits=5), "  [mm]")
						println("			Hpond at end 		= ", ceil(Int, Hpond[Nit]), "  [mm] ")
						println("			∑Runof 				= ", round(∑Runoff_Reduced[end], digits=0), "  [mm] \n" )
		
						println("			Number_of_cells 	        = ", NiZ, "  [-]")

						println("\n			Efficiency 			= ", Efficiency[iOpt_Count], "  [iTer day-1]")
						println("			iNonConverge 			= ", iNonConverge, "  [count]")
						println("			Global_WaterBalance_NormPr 	= ", round(Global_WaterBalance_NormPr[iOpt_Count], digits=8), "  [%]")
						println("			Global_WaterBalance 		= ", 	round(∑∑WaterBalance, digits=8), "  [mm]")
						println("			Average ΔT 			= ",  ΔT_Average[iOpt_Count] , "  [seconds]")
						println("			ΔTmin 				= ",   round(minimum(ΔT[i∑T_CalibrStart:Nit]), digits=0) , "  [seconds]")
						println("			ΔTmax 				= ",  round(maximum(ΔT[i∑T_CalibrStart:Nit]), digits=0) , "  [seconds]")
						println("			ΔT_HyPix 			= ", ceil(Int, ΔRunTimeHypix[iOpt_Count]) , "  [seconds] \n")	


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
						println("			WilmotBest                   = ", round(WilmotBest[iOpt_Count], digits=5))
					end	
					println("		=== === END: summary \n")

					if optionHypix.Table
						tableHypix.TABLE_HYPIX(∑∑ΔSink, ∑WaterBalanceη_Reduced, ∑ΔQ_Bot, CccBest,  Date_Reduced, discret, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, hydroHorizon, iMultistep, iNonConverge_iOpt, N_Layer,Nit_Reduced, NiZ, NseBest, optionHypix, paramHypix, pathOutputHypix, SwcRoots, veg, WilmotBest, WofBest, Z, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced,  ΔPr_Reduced, ΔPrGross_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔRunTimeHypix, ΔSink_Reduced, ΔT_Average, θ_Reduced, θobs_Reduced, θsim_Aver, Ψ_Reduced)
					end
		
					if optionHypix.Ploting
					println("		=== === START: Plotting === ===")
				
						plotHypix.PLOT_HYPIX(∑T_Reduced, clim, Date_Reduced, discret, hydro, hydroHorizon, i∑T_CalibrStart_Day, iMultistep, iScenario, N_iRoot, N_Layer, Nit_Reduced, NiZ, obsTheta, optionHypix, paramHypix, pathOutputHypix, SiteName, veg, Z, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔPrGross_Reduced, ΔQ_Reduced, ΔRootDensity, ΔSink_Reduced, θ_Reduced, θobs_Reduced, θsim_Aver)
		
					println("		=== === END: Plotting === === \n")
					end # if optionHypix.Plotting
			
			println("	=== === === END   ",iMultistep, "  steps ")
			println("	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n")
			end # for loop: iMultistep

		end # for iScenario = 1:N_Scenario
			
	end  # function: HYPIX_START
	# ------------------------------------------------------------------

end  # module hydro
# ............................................................

@time hypixStart.HYPIX_START("LYSIMETERS")
# @time hypixStart.HYPIX_START("TESTCASE")
