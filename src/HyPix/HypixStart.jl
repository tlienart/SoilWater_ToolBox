# =============================================================
#		MODULE: hydro
# =============================================================
module hypixStart

	import ..climate, ..discretization, ..horizonLayer, ..hydroStruct, ..hypixModel, ..hypixOpt, ..interpolate, ..memory, ..ofHypix, ..paths, ..plotHypix, ..reading, ..stats, ..table, ..thetaObs, ..tool, ..vegStruct, ..waterBalance, ..Δtchange, ..θaver
	
	import Statistics, Dates, DelimitedFiles
	export HYPIX_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIX_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYPIX_START(Soilname, option, param,  PathData_Hypix, PathData_SoilWater, SiteName_Hypix, SiteName_Soilwater, Soilwater_OR_Hypix⍰)

		# ===========================================================
		# 					LOOP FOR DIFFERENTY SIMULATIONS
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		# If no optimize
		if !(option.hyPix.Optimisation)
			param.hyPix.iOpt_Start = 1
			param.hyPix.iOpt_End = 1
		end

		for iSim = 1:length(Soilname)
		for iOpt = param.hyPix.iOpt_Start : param.hyPix.iOpt_End
			println("=== START RUNNING Hypix_1D ==== \n")
			println("		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			println("		=== === === START: Looping with time $iOpt steps \n")

			# READING STRUCTURE OF PATH
				path = paths.PATH(iSim, option, PathData_Hypix, PathData_SoilWater, SiteName_Hypix, SiteName_Soilwater, Soilwater_OR_Hypix⍰; Soilname=Soilname)

				println("\n	 ==== ==== ===  $(path.hyPix.IdName_Hypix) 	=== ==== ====\n")

			# COUNT SIMULATIONS
				iOpt_Count = iOpt - param.hyPix.iOpt_Start + 1

			# ACTUAL TIME
				Time_Start = Dates.now()

			# DISCRETISATION ~~~~~
				# Read discretisaation
					Layer, N_iHorizon, N_iZ, Z, θ_Ini = reading.hyPix.DISCRETIZATION(path.hyPix)
				# Process discretisation ~~~~~
					discret = discretization.DISCRETIZATION(N_iZ, Z)

			# CLIMATE DATA  ~~~~~
				# Read climate
					clim = reading.hyPix.CLIMATE(option, param, path.hyPix)
				# Process climate
					∑Pet_Climate, ∑Pr_Climate, ∑T_Climate, N_∑T_Climate, Temp = climate.CLIMATE(clim, option)

			# OBSERVED θ  ~~~~~
			if option.hyPix.θobs
				# Read observed θ
					obsTheta = reading.hyPix.TIME_SERIES(option, param, path.hyPix)
				
				# Process observed θ
					obsTheta = thetaObs.ΘOBS(obsTheta, clim, discret, Z)
			end #  option.hyPix.θobs 

			# MEMORY
				if iOpt_Count == 1
					N_∑T_Plot = param.hyPix.iOpt_End - param.hyPix.iOpt_Start + 1

					global Efficiency                 = fill(0.0::Float64, N_∑T_Plot)
					global Global_WaterBalance        = fill(0.0::Float64, N_∑T_Plot)
					global Global_WaterBalance_NormPr = fill(0.0::Float64, N_∑T_Plot)
					global RmseBest                   = fill(0.0::Float64, N_∑T_Plot)
					global SwcRoots                   = fill(0.0::Float64, N_∑T_Plot)
					global WofBest                    = fill(0.0::Float64, N_∑T_Plot)
					global ΔRunTimeHypix              = fill(0.0::Float64, N_∑T_Plot)
					global ΔT_Average                 = fill(0.0::Float64, N_∑T_Plot)
					global ∑ΔQ_Bot                    = fill(0.0::Float64, N_∑T_Plot)
					global ∑∑ΔSink                    = fill(0.0::Float64, N_∑T_Plot)
				end
				∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pr, ∑T, CropCoeficientᵀ, iNonConverge_iOpt, Laiᵀ, N_Memory, Q, Residual, ΔEvaporation, ΔHpond, ΔΨmax, ΔPet, ΔPr, ΔSink, ΔT, θ, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest = memory.MEMORY(clim, iOpt_Count, N_∑T_Climate, N_iZ, obsTheta, param)

			# INITIALIZING THE STRUCTURE
				# Initializing hydroHorizon structure
					hydroHorizon = hydroStruct.HYDROSTRUCT(option.hyPix, N_iHorizon)

				# Initializing hydraulic param into structure 
					hydro = hydroStruct.HYDROSTRUCT(option.hyPix, N_iZ)

				# Initialiozing  vegetation parameters into veg structure
					veg = vegStruct.VEGSTRUCT()

				# Optimisation
					hydroHorizon_best = hydroStruct.HYDROSTRUCT(option.hyPix, N_iHorizon)
					hydro_best        = hydroStruct.HYDROSTRUCT(option.hyPix, N_iZ)
					veg_best          = vegStruct.VEGSTRUCT()
				
			# ========================================================
	 
			# OBTAINING HYDRAULIC AND VEGETATION PARAMETERS
			if option.hyPix.Optimisation
				hydro, hydroHorizon, optim, veg = reading.hyPix.HYPIX_PARAM(Layer, hydro, hydroHorizon, iOpt, N_iZ, option, param, path.hyPix.HyPixParamOpt, veg)
			else
			# Reading veg parameters
				veg, ~ = reading.READ_STRUCT(veg, path.hyPix.HyPix_VegParam)
	
			# Hydraulic parameters
				hydroHorizon, ~ = reading.READ_STRUCT(hydroHorizon, path.hyPix.HyPix_HydroParam)

				hydro        = horizonLayer.HYDROHORIZON_2_HYDRO(hydroHorizon, Layer, N_iZ, option)

			# options of optim		
				Flag_Opt = false
				NparamOpt = 0

				optim = ( NparamOpt=NparamOpt, Flag_Opt=Flag_Opt)		
			end # option.hyPix.Optimisation

			# SINK TERM 
            Laiᵀ_η            = reading.hyPix.LOOKUPTABLE_LAI(clim, option, path.hyPix, veg)
            CropCoeficientᵀ_η = reading.hyPix.LOOKUPTABLE_CROPCOEFICIENT(clim, option, path.hyPix, veg)

			if optim.Flag_Opt
				hydro, hydro_best, hydroHorizon, hydroHorizon_best, veg, veg_best, WofBest = hypixOpt.HYPIXOPTIMISATION_START(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, hydro, hydro_best, hydroHorizon, hydroHorizon_best, iOpt_Count, Laiᵀ, Laiᵀ_η, Layer, N_∑T_Climate, N_iHorizon, N_iZ, obsTheta, optim, option, param, Q, Residual, veg, veg_best, WofBest, Z, ΔEvaporation, ΔHpond, ΔPet, ΔPr, ΔSink, ΔT, ΔΨmax, θ, θ_Ini, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest)
			end
			
			# if Flag_Opt then it will rerun with the optimal parameters
			∑Pet, ∑Pr, ∑T, ∑T_Climate, clim, discret, iNonConverge, Iter_CountTotal, N_iRoot, N_iT, N_iZ, Q, veg, ΔEvaporation, ΔHpond, ΔRootDensity, ΔT, θ, Ψ = hypixModel.HYPIX(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, hydro, Laiᵀ, Laiᵀ_η, N_∑T_Climate, N_iZ, option, param, Q, Residual, veg, Z, ΔEvaporation, ΔHpond, ΔPet, ΔPr, ΔSink, ΔT, ΔΨmax, θ, θ_Ini, Ψ, Ψ_Max, Ψ_Min, Ψbest)

			

			# WATER BALANCE
				# Computed after the warmup period
				∑∑WaterBalance, ∑WaterBalance_η, ∑ΔSink, i∑T_CalibrStart, ΔStorage = waterBalance.WATERBALANCE(∑T, obsTheta, discret, hydro, N_iRoot, N_iT, N_iZ, Q, ΔSink, ΔT, θ, Ψ)


			# SUMMARY HOW GOOD THE SIMULATION
				# Computed climate day after the warmup period
				i∑T_CalibrStart_Day = 1::Int64 
				while ∑T_Climate[i∑T_CalibrStart_Day] < obsTheta.∑T[1]
					i∑T_CalibrStart_Day += 1
				end
				i∑T_CalibrStart_Day += 1
				∑Pr_Clim = sum(clim.Pr[i∑T_CalibrStart_Day:clim.N_Climate]) 
				∑Pet_Net = sum(clim.Pet[i∑T_CalibrStart_Day:clim.N_Climate])

				# Soil water content for the rootzone at the end of simulation
				@fastmath @inbounds @simd for iZ=1:N_iRoot
					SwcRoots[iOpt_Count] += θ[N_iT,iZ] * discret.ΔZ[iZ]
				end

				# Timing 
					Time_End = Dates.now()
					ΔRunTimeHypix[iOpt_Count] = Dates.value(Time_End - Time_Start) / 1000

				# Convergence rate
					iNonConverge_iOpt[iOpt_Count]          = iNonConverge

					Efficiency[iOpt_Count]                 = ceil(Int, 86400 * Iter_CountTotal / ∑T[N_iT])
					
					ΔT_Average[iOpt_Count]                 = ceil(Int, Statistics.mean(ΔT[i∑T_CalibrStart:N_iT]))
					
					Global_WaterBalance[iOpt_Count]        = ∑∑WaterBalance
					
					Global_WaterBalance_NormPr[iOpt_Count] = 100.0 * ∑WaterBalance_η[N_iT]
					
					∑∑ΔSink[iOpt_Count]                    = ∑ΔSink[N_iT]
				
				# Ground water recharge
					∑ΔQ_Bot[iOpt_Count] = 0.0
					@fastmath @inbounds @simd for iT=i∑T_CalibrStart:N_iT
						∑ΔQ_Bot[iOpt_Count] += ΔT[iT] * Q[iT, N_iZ+1]
					end
           
				println("		=== ===START: summary  $iOpt steps ...")

					println("			∑Pr 			= ", ceil(Int, ∑Pr_Clim), "  [mm]")
					println("			∑Pr_Soil 		= ", ceil(Int, ∑Pr[N_iT] - ∑Pr[i∑T_CalibrStart]),  "  [mm]")
					println("			∑Pr_Intercepted/∑Pr 	= ", ceil(Int, 100. * (∑Pr_Clim - (∑Pr[N_iT]-∑Pr[i∑T_CalibrStart])) / ∑Pr_Clim),  "  [%]")
					println("			∑Pet_Net 		= ", ceil(Int, ∑Pet_Net), "  [mm]")
					println("			∑Pet 			= ", ceil(Int, ∑Pet[N_iT]- ∑Pet[i∑T_CalibrStart]), "  [mm]")
					println("			∑ΔSink/∑Pet_Net 	= ", ceil(Int, 100.0 * ∑ΔSink[N_iT] / ∑Pet_Net), "  [%] \n")
					
					println("			∑SoilWaterContentRootEnd = ", ceil(Int, SwcRoots[iOpt_Count]), "  [mm]")
					println("			∑ΔSink 			= ", -ceil(Int, ∑∑ΔSink[iOpt_Count]), "  [mm]")
					println("			∑Infilt_Bot 		= ", -ceil( Int, ∑ΔQ_Bot[iOpt_Count] ), "  [mm]")
					println("			ΔHpond at end 		= ", ceil(Int, ΔHpond[N_iT]), "  [mm] \n")

					println("			iNonConverge 			= ", iNonConverge_iOpt[iOpt_Count], "  [count]")
					println("			Global_WaterBalance_NormPr 	= ", round(Global_WaterBalance_NormPr[iOpt_Count], digits=2), "  [%]")
					println("			Global_WaterBalance 		= ", 	round(Global_WaterBalance[iOpt_Count], digits=2), "  [mm]")
					println("			Average ΔT 			= ",  ΔT_Average[iOpt_Count] , "  [seconds]")
					println("			ΔTmin 				= ",   round(minimum(ΔT[i∑T_CalibrStart:N_iT]), digits=0) , "  [seconds]")
					println("			ΔTmax 				= ",  round(maximum(ΔT[i∑T_CalibrStart:N_iT]), digits=0) , "  [seconds]")
					println("			ΔT_HyPix 			= ", ceil(Int, ΔRunTimeHypix[iOpt_Count]) , "  [seconds]")			
					println("			Efficiency 			= ", Efficiency[iOpt_Count], "  [iTer day-1], \n")


				∑T_Plot, ∑T_Date_Plot, ∑WaterBalance_η_Plot, Date_Plot, Flag_Plot_Pond, N_∑T_Plot, ΔEvaporation_Plot, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, ΔT_Plot, θ_Plot, θobs_Plot, Ψ_Plot = Δtchange.CHANGE_OUTPUT_ΔT(∑Pet[1:N_iT], ∑Pr[1:N_iT], ∑T[1:N_iT], ∑WaterBalance_η[1:N_iT], ∑ΔSink[1:N_iT], obsTheta, clim, N_iT, N_iZ, param, Q[1:N_iT,1:N_iZ+1], veg, ΔEvaporation[1:N_iT], ΔHpond[1:N_iT], ΔT[1:N_iT], θ[1:N_iT,1:N_iZ], Ψ[1:N_iT,1:N_iZ], ∑T_Climate, path.hyPix)

			# Computing average simulated θ to comapre it with average observed θ
			if option.hyPix.θobs_Average && option.hyPix.θobs	
				θsim_Aver = θaver.θAVER(discret; Z=Z, θ_Plot=θ_Plot, N_iZ=N_iZ, N_∑T_Plot=N_∑T_Plot, Zaver=400.0)

				RmseBest[iOpt_Count] = stats.NSE(θobs_Plot[1:N_∑T_Plot], θsim_Aver[1:N_∑T_Plot])
				
			elseif !(option.hyPix.θobs) && option.hyPix.θobs
				RmseBest[iOpt_Count] = ofHypix.θof.RMSE_θ(∑T, obsTheta, N_iT, N_iZ, θ, θSim)
				θsim_Aver = []
			end

			if option.hyPix.θobs
				println("			Nse 			= ", round(RmseBest[iOpt_Count], digits=5), "  [mm3 mm-3]")
			end
			println("		=== === END: summary \n")


			if option.hyPix.Table
			println("		=== === START: Table === ===")
	
				# Writing values of hydraulic parameters
				table.hyPix.HYDRO(hydroHorizon, iOpt, N_iHorizon, path.hyPix)

				# Writing values of veg parameters
				table.hyPix.VEG(veg, iOpt, path.hyPix)

				table.hyPix.PERFORMANCE(∑∑ΔSink, ∑ΔQ_Bot, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iOpt, iOpt, RmseBest, SwcRoots, WofBest, ΔRunTimeHypix, ΔT_Average, Soilname, path.hyPix)		

				if option.hyPix.Table_Discretization
					table.hyPix.DISCRETIZATION(discret, N_iZ, Z[1:N_iZ], path.hyPix)
				end
				if  option.hyPix.Table_TimeSeries
					table.hyPix.TIME_SERIES(∑T[1:N_iT], ΔT[1:N_iT], ∑Pr[1:N_iT], ΔPr[1:N_iT], ΔHpond[1:N_iT], ΔT[1:N_iT].*Q[1:N_iT,1], ∑WaterBalance_η[1:N_iT], iOpt, path.hyPix)
				end
				if option.hyPix.Table_TimeSeriesDaily
					table.hyPix.TIME_SERIES_DAILY(∑T_Plot[1:N_∑T_Plot], ∑WaterBalance_η_Plot[1:N_∑T_Plot], Date_Plot[1:N_∑T_Plot], iOpt, N_∑T_Plot, ΔEvaporation_Plot[1:N_∑T_Plot], ΔFlux_Plot[1:N_∑T_Plot, N_iZ+1], ΔPet_Plot[1:N_∑T_Plot], ΔPond_Plot[1:N_∑T_Plot], ΔPr_Plot[1:N_∑T_Plot], ΔSink_Plot[1:N_∑T_Plot], path.hyPix)
				end
				if option.hyPix.Table_θ
					table.hyPix.θ(∑T[1:N_iT], θ[1:N_iT,1:N_iZ], discret.Znode[1:N_iZ], iOpt, path.hyPix)
				end
				if option.hyPix.Table_Ψ
					table.hyPix.Ψ(∑T[1:N_iT], Ψ[1:N_iT,1:N_iZ], discret.Znode[1:N_iZ], iOpt, path.hyPix)
				end
				if option.hyPix.Table_Q
					table.hyPix.Q(∑T[1:N_iT], Q[1:N_iT,1:N_iZ+1], Z[N_iZ], discret.Znode[1:N_iZ], iOpt, path.hyPix)
				end
				if option.hyPix.Tabule_θΨ
					table.hyPix.θΨ(hydroHorizon, iOpt, N_iHorizon, option.hyPix, param, path.hyPix)
					table.hyPix.KΨ(hydroHorizon, iOpt, N_iHorizon, option.hyPix, param, path.hyPix)
				end
				if option.hyPix.Table_Climate
					table.hyPix.DAILY_CLIMATE(∑T_Climate, clim, iOpt, path.hyPix)
				end
				if option.hyPix.θobs_Average && option.hyPix.θobs
					table.hyPix.θAVERAGE(Date_Plot[1:N_∑T_Plot], iOpt, θobs_Plot[1:N_∑T_Plot], θsim_Aver[1:N_∑T_Plot], path.hyPix)
				end
			println("		=== === END: Table === === \n")
			end  # if option.hyPix.Table
				
			if option.other.Ploting
			println("		=== === START: Plotting === ===")

				# if option.hyPix.Plot_Other
				# 	# plotOther.ΨMINΨMAX(hydro, path.hyPix)
				# 	plotOther.WOF_STEPS(path.hyPix)
				# 	# plotOther.SE_Ψ_CONSTRAINED(hydro, path.hyPix)
				# 	# plotOther.PLOT_σ_2_θr(hydro, path.hyPix)
				# 	# plotOther.PLOT_θΨ_Δθ(hydro, path.hyPix)
				# end # option.hyPix.Plot_Other
				
				if option.hyPix.Plot_Hypix
					plotHypix.plots.TIMESERIES(∑T_Date_Plot, ∑T_Plot, obsTheta, discret, Flag_Plot_Pond, iOpt, N_∑T_Plot, N_iZ, option, param, ΔEvaporation_Plot, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, θ_Plot, θobs_Plot, clim, i∑T_CalibrStart_Day, θsim_Aver, path.hyPix)
					# plotHypix.TIME_SERIES(∑T_Plot, ∑WaterBalance_η_Plot, obsTheta, discret, Flag_Plot_Pond, iOpt, N_∑T_Plot, N_iZ, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, ΔT_Plot, θ_Plot, θobs_Plot, Ψ_Plot, path.hyPix)
				end
				if option.hyPix.Plot_θΨK
					plotHypix.θΨK(hydroHorizon, N_iHorizon, iOpt, path.hyPix)
				end
				if option.hyPix.Plot_Vegetation && option.hyPix.RootWaterUptake
					plotHypix.VEG_FUNCTIONS(discret, iOpt, N_iRoot, veg, Z, ΔRootDensity, path.hyPix)
				end
				if option.hyPix.Plot_Interception
					plotHypix.plots.RAINFALL_INTERCEPTION(clim, i∑T_CalibrStart_Day, iOpt, path.hyPix)
				end
				if  option.hyPix.Plot_Sorptivity
					plotHypix.plots.PLOT_SORPTIVITY(hydro, iOpt, option, optionₘ, path.hyPix)
				end
			println("		=== === END: Plotting === === \n")
			end # if option.hyPix.Plotting
	
			println("	=== === === END  : Looping with time ")
			println("	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n")

	# Garbage collect
	GC.gc()
	end # for loop: iOpt
	end #iSim = 1:length(Soilname)

	println("\n ==== END RUNNING Hypix_1D =====")

	return
	end # function HYPIX_START
	
end  # module hydro
# ............................................................