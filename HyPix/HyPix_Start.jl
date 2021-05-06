# =============================================================
#		MODULE: hydro
# =============================================================
module hypix

	import ..hydroStruct, ..hypixModel, ..hypixOpt, ..interpolate, ..memory, ..ofHypix, ..option, ..param, ..path, ..plotHypix, ..plotOther, ..priorprocess, ..readHypix, ..signature, ..tableHypix, ..tool, ..vegStruct, ..waterBalance, ..zobs, ..Δtchange, ..discretization
	using Statistics, Dates
	export HYPIX_START		

	function HYPIX_START()
		println("=== START RUNNING Hypix_1D ==== \n")

		Time_Start = Dates.now()

		println("\n	==== ==== ===  $(path.SiteName_Hypix) 	=== ==== ====\n")

		Horizon, N_iHorizon, N_iZ, Z, θ_Ini = readHypix.DISCRETIZATION()
	
		# Discretisation ~~~~~
			discret = discretization.DISCRETIZATION(N_iZ, Z)

		# Climate data  ~~~~~
			clim = readHypix.CLIMATE()

			∑Pet_Climate, ∑Pr_Climate, ∑T_Climate, N_∑T_Climate, Temp = priorprocess.CLIMATE(clim)

		# calibr data  ~~~~~
			if option.hypix.calibr
				calibr = readHypix.TIME_SERIES()
				
				# Deriving calibr.∑T and Cells for which calibr are measured
				calibr = zobs.ZOBS(calibr, clim, discret, Z)
			end #  option.hypix.calibr

			∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑∑ΔSink, ∑Pet, ∑Pr, ∑T, ∑ΔQ_Bot, CropCoeficientᵀ, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iSim, Laiᵀ, N_Memory, Q, Residual, RmseBest, SwcRoots, WofBest, ΔEvaporation, ΔHpond, ΔΨmax, ΔPet, ΔPr, ΔRunTimeHypix, ΔSink, ΔT, ΔT_Average, θ, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest = memory.MEMORY(calibr, clim, N_∑T_Climate, N_iZ)
		
		# INITIALIZING THE STRUCTURE BEFORE OPTIMISATION
			# Initializing the hydraulic parameters into hydroHorizon structure
				hydroHorizon = hydroStruct.HYDROSTRUCT(N_iHorizon)

			# Initializing hydraulic param into hydro structure 
				hydro = hydroStruct.HYDROSTRUCT(N_iZ)

			# Initialiozing the vegetation parameters into veg structure
				veg = vegStruct.VEGSTRUCT()

			# For optimisation
				hydroHorizon_best = hydroStruct.HYDROSTRUCT(N_iHorizon)
				hydro_best        = hydroStruct.HYDROSTRUCT(N_iZ)
				veg_best          = vegStruct.VEGSTRUCT()
	 
		# =============================================
		# Loop for the different simulations loop
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for iSim = param.hypix.iSim_Start : param.hypix.iSim_End
		println("	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		println("	=== === === START: Looping with time $iSim steps \n")

			iSim_Count = iSim - param.hypix.iSim_Start + 1

			# Reading hydraulic and veg param ~~~~~
				hydro, hydroHorizon, optim, veg = readHypix.HYPIX_PARAM(Horizon, hydro, hydroHorizon, iSim, N_iZ, veg)

			# SINK TERM 
            Laiᵀ_η            = readHypix.LOOKUPTABLE_LAI(clim, veg)
            CropCoeficientᵀ_η = readHypix.LOOKUPTABLE_CROPCOEFICIENT(clim, veg)

			if optim.Flag_Opt
				hydro, hydro_best, hydroHorizon, hydroHorizon_best, veg, veg_best, WofBest = hypixOpt.HYPIXOPTIMISATION_START(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, calibr, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Horizon, hydro, hydro_best, hydroHorizon, hydroHorizon_best, iSim_Count, Laiᵀ, Laiᵀ_η, N_∑T_Climate, N_iHorizon, N_iZ, optim, Q, Residual, veg, veg_best, WofBest, Z, ΔEvaporation, ΔHpond, ΔΨmax, ΔPet, ΔPr, ΔSink, ΔT, θ, θ_Ini, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest)
			end
			
			# if Flag_Opt then it will rerun with the optimal parameters
			∑Pet, ∑Pr, ∑T, ∑T_Climate, clim, discret, iNonConverge, Iter_CountTotal, N_iRoot, N_iT, N_iZ, Q, veg, ΔEvaporation, ΔHpond, ΔRootDensity, ΔT, θ, Ψ = hypixModel.HYPIX(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, hydro, Laiᵀ, Laiᵀ_η, N_∑T_Climate, N_iZ, Q, Residual, veg, Z, ΔEvaporation, ΔHpond, ΔΨmax, ΔPet, ΔPr, ΔSink, ΔT, θ, θ_Ini, Ψ, Ψ_Max, Ψ_Min, Ψbest)

			# WATER BALANCE
			# Computed after the warmup period
			∑∑WaterBalance, ∑WaterBalance_η, ∑ΔSink, i∑T_CalibrStart, ΔStorage = waterBalance.WATERBALANCE(∑T, calibr, discret, hydro, N_iRoot, N_iT, N_iZ, Q, ΔSink, ΔT, θ, Ψ)

			# SIGNATURE
			if option.hypix.Signature_Run
				Signature_Deficit_Obs, Signature_Max_Obs, Signature_Saturated_Obs, Signature_Senescence_Obs, Signature_Deficit_Sim, Signature_Max_Sim, Signature_Saturated_Sim, Signature_Senescence_Sim = signature.SIGNATURE(∑T[1:N_iT], calibr, hydroHorizon, N_iRoot, N_iT, N_iZ, veg, ΔRootDensity, Ψ[1:N_iT,1:N_iZ])
			end 

			# SUMMARY HOW GOOD THE SIMULATION
				# Computed climate day after the warmup period
				i∑T_CalibrStart_Day = 1::Int64 
				while ∑T_Climate[i∑T_CalibrStart_Day] < calibr.∑T[1]
					i∑T_CalibrStart_Day += 1
				end
				i∑T_CalibrStart_Day += 1
				∑Pr_Clim = sum(clim.Pr[i∑T_CalibrStart_Day:clim.N_Climate]) 
				∑Pet_Net = sum(clim.Pet[i∑T_CalibrStart_Day:clim.N_Climate])
				
				# In optim.Flag_Opt it is already computed
				if option.hypix.calibr
					RmseBest[iSim_Count] = ofHypix.θof.RMSE_θ(∑T, calibr, N_iT, N_iZ, θ, θSim)
				end

				# Soil water content for the rootzone at the end of simulation
				@fastmath @inbounds @simd for iZ=1:N_iRoot
					SwcRoots[iSim_Count] += θ[N_iT,iZ] * discret.ΔZ[iZ]
				end

				# Timing 
					Time_End = Dates.now()
					ΔRunTimeHypix[iSim_Count] = Dates.value(Time_End - Time_Start) / 1000

				# Convergence rate
					iNonConverge_iSim[iSim_Count]          = iNonConverge

					Efficiency[iSim_Count]                 = ceil(Int, 86400 * Iter_CountTotal / ∑T[N_iT])
					
					ΔT_Average[iSim_Count]                 = ceil(Int, Statistics.mean(ΔT[i∑T_CalibrStart:N_iT]))
					
					Global_WaterBalance[iSim_Count]        = ∑∑WaterBalance
					
					Global_WaterBalance_NormPr[iSim_Count] = 100.0 * ∑WaterBalance_η[N_iT]
					
					∑∑ΔSink[iSim_Count]                    = ∑ΔSink[N_iT]
				
				# Ground water recharge
					∑ΔQ_Bot[iSim_Count] = 0.0
					@fastmath @inbounds @simd for iT=i∑T_CalibrStart:N_iT
						∑ΔQ_Bot[iSim_Count] += ΔT[iT] * Q[iT, N_iZ+1]
					end
           
				println("		=== ===START: summary  $iSim steps ...")

					println("			∑Pr 			= ", ceil(Int, ∑Pr_Clim), "  [mm]")
					println("			∑Pr_Soil 		= ", ceil(Int, ∑Pr[N_iT] - ∑Pr[i∑T_CalibrStart]),  "  [mm]")
					println("			∑Pr_Intercepted/∑Pr 	= ", ceil(Int, 100. * (∑Pr_Clim - (∑Pr[N_iT]-∑Pr[i∑T_CalibrStart])) / ∑Pr_Clim),  "  [%]")
					println("			∑Pet_Net 		= ", ceil(Int, ∑Pet_Net), "  [mm]")
					println("			∑Pet 			= ", ceil(Int, ∑Pet[N_iT]- ∑Pet[i∑T_CalibrStart]), "  [mm]")
					println("			∑ΔSink/∑Pet_Net 	= ", ceil(Int, 100.0 * ∑ΔSink[N_iT] / ∑Pet_Net), "  [%] \n")
					
					println("			∑SoilWaterContentRootEnd = ", ceil(Int, SwcRoots[iSim_Count]), "  [mm]")
					println("			∑ΔSink 			= ", -ceil(Int, ∑∑ΔSink[iSim_Count]), "  [mm]")
					println("			∑Infilt_Bot 		= ", -ceil( Int, ∑ΔQ_Bot[iSim_Count] ), "  [mm]")
					println("			ΔHpond at end 		= ", ceil(Int, ΔHpond[N_iT]), "  [mm] \n")

					println("			iNonConverge 			= ", iNonConverge_iSim[iSim_Count], "  [count]")
					println("			Global_WaterBalance_NormPr 	= ", round(Global_WaterBalance_NormPr[iSim_Count], digits=2), "  [%]")
					println("			Global_WaterBalance 		= ", 	round(Global_WaterBalance[iSim_Count], digits=2), "  [mm]")
					println("			Average ΔT 			= ",  ΔT_Average[iSim_Count] , "  [seconds]")
					println("			ΔTmin 				= ",   round(minimum(ΔT[i∑T_CalibrStart:N_iT]), digits=0) , "  [seconds]")
					println("			ΔTmax 				= ",  round(maximum(ΔT[i∑T_CalibrStart:N_iT]), digits=0) , "  [seconds]")
					println("			ΔT_HyPix 			= ", ceil(Int, ΔRunTimeHypix[iSim_Count]) , "  [seconds]")			
					println("			Efficiency 			= ", Efficiency[iSim_Count], "  [iTer day-1], \n")

					if option.hypix.calibr
						println("			RmseBest 			= ", round(RmseBest[iSim_Count], digits=5), "  [mm3 mm-3]")
					end
				println("		=== === END: summary \n")


				println("		=== === START: increasing time step === ===")

				∑T_Plot, ∑T_Date_Plot, ∑WaterBalance_η_Plot, Date_Plot, Flag_Plot_Pond, N_∑T_Plot, ΔEvaporation_Plot, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, ΔT_Plot, θ_Plot, θobs_Plot, Ψ_Plot = Δtchange.CHANGE_OUTPUT_ΔT(∑Pet[1:N_iT], ∑Pr[1:N_iT], ∑T[1:N_iT], ∑WaterBalance_η[1:N_iT], ∑ΔSink[1:N_iT], calibr, clim, N_iT, N_iZ, Q[1:N_iT,1:N_iZ+1], veg, ΔEvaporation[1:N_iT], ΔHpond[1:N_iT], ΔT[1:N_iT], θ[1:N_iT,1:N_iZ], Ψ[1:N_iT,1:N_iZ], ∑T_Climate)
			
				println("		=== === END: increasing time step === ===")

				
			if option.Plot
			println("		=== === START: Plotting === ===")

				if option.hypix.Plot_Other
					# plotOther.ΨMINΨMAX(hydro)
					plotOther.WOF_STEPS()
					# plotOther.SE_Ψ_CONSTRAINED(hydro)
					# plotOther.PLOT_σ_2_θr(hydro)
					# plotOther.PLOT_θΨ_Δθ(hydro)
				end # option.hypix.Plot_Other
				
				if option.hypix.Plot_Hypix

					plotHypix.plots.TIMESERIES(∑T_Date_Plot, ∑T_Plot, calibr, discret, Flag_Plot_Pond, iSim, N_∑T_Plot, N_iZ, ΔEvaporation_Plot, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, θ_Plot, θobs_Plot, clim, i∑T_CalibrStart_Day)
					
					# plotHypix.TIME_SERIES(∑T_Plot, ∑WaterBalance_η_Plot, calibr, discret, Flag_Plot_Pond, iSim, N_∑T_Plot, N_iZ, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, ΔT_Plot, θ_Plot, θobs_Plot, Ψ_Plot)
				end
				if option.hypix.Plot_θΨK
					plotHypix.θΨK(hydroHorizon, N_iHorizon, iSim)
				end
				if option.hypix.Plot_Vegetation && option.hypix.RootWaterUptake
					plotHypix.VEG_FUNCTIONS(discret, iSim, N_iRoot, veg, Z, ΔRootDensity)
				end
				if option.hypix.Plot_Interception
					plotHypix.plots.RAINFALL_INTERCEPTION(clim, i∑T_CalibrStart_Day, iSim)
				end
				if  option.hypix.Plot_Sorptivity
					plotHypix.plots.PLOT_SORPTIVITY(iSim, hydro)
				end
			println("		=== === END: Plotting === === \n")
			end # if option.hypix.Plotting

			
			if option.hypix.Table
			println("		=== === START: Table === ===")

				# Writing values of hydraulic parameters
				tableHypix.HYDRO(hydroHorizon, iSim, N_iHorizon)

				# Writing values of veg parameters
				tableHypix.VEG(veg, iSim)

				tableHypix.PERFORMANCE(∑∑ΔSink, ∑ΔQ_Bot, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iSim, iSim, RmseBest, SwcRoots, WofBest, ΔRunTimeHypix, ΔT_Average)		

				if option.hypix.Table_Discretization
					tableHypix.DISCRETIZATION(discret, N_iZ, Z[1:N_iZ])
				end
				if  option.hypix.Table_TimeSeries
					tableHypix.TIME_SERIES(∑T[1:N_iT], ΔT[1:N_iT], ∑Pr[1:N_iT], ΔPr[1:N_iT], ΔHpond[1:N_iT], ΔT[1:N_iT].*Q[1:N_iT,1], ∑WaterBalance_η[1:N_iT], iSim)
				end
				if option.hypix.Table_TimeSeriesDaily
					tableHypix.TIME_SERIES_DAILY(∑T_Plot[1:N_∑T_Plot], ∑WaterBalance_η_Plot[1:N_∑T_Plot], Date_Plot[1:N_∑T_Plot], iSim, N_∑T_Plot, ΔEvaporation_Plot[1:N_∑T_Plot], ΔFlux_Plot[1:N_∑T_Plot, N_iZ+1], ΔPet_Plot[1:N_∑T_Plot], ΔPond_Plot[1:N_∑T_Plot], ΔPr_Plot[1:N_∑T_Plot], ΔSink_Plot[1:N_∑T_Plot])
				end
				if option.hypix.Table_θ
					tableHypix.θ(∑T[1:N_iT], θ[1:N_iT,1:N_iZ], discret.Znode[1:N_iZ], iSim)
				end
				if option.hypix.Table_Ψ
					tableHypix.Ψ(∑T[1:N_iT], Ψ[1:N_iT,1:N_iZ], discret.Znode[1:N_iZ], iSim)
				end
				if option.hypix.Table_Q
					tableHypix.Q(∑T[1:N_iT], Q[1:N_iT,1:N_iZ+1], Z[N_iZ], discret.Znode[1:N_iZ], iSim)
				end
				if option.hypix.Signature_Run
					tableHypix.SIGNATURE(iSim, Signature_Deficit_Obs, Signature_Max_Obs, Signature_Saturated_Obs, Signature_Senescence_Obs, Signature_Deficit_Sim, Signature_Max_Sim, Signature_Saturated_Sim, Signature_Senescence_Sim)
				end
				if option.hypix.Tabule_θΨ
					tableHypix.θΨ(hydroHorizon, iSim, N_iHorizon)
					tableHypix.KΨ(hydroHorizon, iSim, N_iHorizon)
				end
				if option.hypix.Table_Climate
					tableHypix.DAILY_CLIMATE(∑T_Climate, clim, iSim)
				end
			println("		=== === END: Table === === \n")
			end  # if option.hypix.Table

	println("	=== === === END  : Looping with time ")
	println("	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n")
	end # for loop: iSim

	println("\n ==== END RUNNING Hypix_1D =====")
	return
	end # function HYPIX_START
	
end  # module hydro
# ............................................................