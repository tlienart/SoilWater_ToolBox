# =============================================================
#		MODULE: hydro
	# =============================================================

	include("Interpolate.jl")
	include("Opt/ThetaObs.jl")
	include("θini.jl")
	include("Opt/OfHypix.jl")
	include("Interception.jl")
	include("Flux.jl")
	include("Ponding.jl")
	include("Residual.jl")
	include("ΔΔtchange.jl")
	include("TimeStep.jl")
	include("Richard.jl")
	include("WaterBalance.jl")
	include("Evaporation.jl")
	include("RootWaterUptake.jl")
	include("CheckError.jl")
	include("Pet.jl")
	include("Other/θaver.jl")
	include("Memory.jl")
	include("Climate.jl")
	include("PlotHypix.jl")
	include("HypixModel.jl")
	include("Opt/HypixOpt.jl")
	include("Other/PlotOther.jl")


module hypixStart
	import ..climate, ..cst, ..discretization, ..horizonLayer, ..hydroStruct, ..hypixModel, ..hypixOpt, ..interpolate, ..memory, ..ofHypix, ..paths, ..plotHypix, ..reading, ..stats, ..table, ..thetaObs, ..tool, ..vegStruct, ..waterBalance, ..Δtchange, ..θaver, ..plotOther
	import Statistics: mean
	import Dates: now, value
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
		
		# LOOPING FOR EVERY SOIL	
		for iSim = 1:length(Soilname)

			Nit_Reduced = param.hyPix.iOpt_End - param.hyPix.iOpt_Start + 1

			Efficiency                 = fill(0.0::Float64, Nit_Reduced)
			Global_WaterBalance        = fill(0.0::Float64, Nit_Reduced)
			Global_WaterBalance_NormPr = fill(0.0::Float64, Nit_Reduced)
			RmseBest                   = fill(0.0::Float64, Nit_Reduced)
			SwcRoots                   = fill(0.0::Float64, Nit_Reduced)
			WofBest                    = fill(0.0::Float64, Nit_Reduced)
			ΔRunTimeHypix              = fill(0.0::Float64, Nit_Reduced)
			ΔT_Average                 = fill(0.0::Float64, Nit_Reduced)
			∑ΔQ_Bot                    = fill(0.0::Float64, Nit_Reduced)
			∑∑ΔSink                    = fill(0.0::Float64, Nit_Reduced)

		for iOpt = param.hyPix.iOpt_Start:param.hyPix.iOpt_End
			println("=== START RUNNING Hypix_1D ==== \n")
			println("		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			println("		=== === === START: Looping with time $iOpt steps \n")

			# READING STRUCTURE OF PATH
				path = paths.PATH(iSim, option, PathData_Hypix, PathData_SoilWater, SiteName_Hypix, SiteName_Soilwater, Soilwater_OR_Hypix⍰; Soilname=Soilname)

				println("\n	 ==== ==== ===  $(path.hyPix.IdName_Hypix) 	=== ==== ====\n")

			# COUNT SIMULATIONS
				iOpt_Count = iOpt - param.hyPix.iOpt_Start + 1

			# ACTUAL TIME
				Time_Start = now()

			# DATES OF SIMULATION
				param = reading.hyPix.DATES(param, path.hyPix)

			# DISCRETISATION ~~~~~
			# Create Discretisation.csv from SoilLayer.csv
			if option.hyPix.Discretisation_File_Auto⍰ == "Auto"
				# Read SoilLayer, could be either θini, Ψini
				Flag_θΨini, Layer, N_Layer, ~, Zlayer, θini, Ψini = reading.hyPix.DISCRETIZATION(path.hyPix.DiscretizationAuto)

				# Performing y6auto discretisation			
					Layer, Z, θΨini_Cell = discretization.DISCRETIZATION_AUTO(option, param; Flag_θΨini=Flag_θΨini, N_Layer=N_Layer, Zlayer=Zlayer, θini=θini, Ψini=Ψini)
		
					table.hyPix.DISCRETIZATION_AUTO(Flag_θΨini, Layer, path.hyPix, Z, θΨini_Cell)

					@info "***			Created new Discretisation.csv file ***" 
			end # if option.hyPix.Discretisation_File_Auto⍰ == "Auto" 
		
			# Read discretisation
				Flag_θΨini, Layer, N_Layer, NiZ, Z, θini, Ψini = reading.hyPix.DISCRETIZATION(path.hyPix.Discretization)

			# Process discretisation of the soil profile ~~~~~
				discret = discretization.DISCRETIZATION(NiZ, Z)

			# CLIMATE DATA  ~~~~~
				# Read climate
					clim = reading.hyPix.CLIMATE(option, param, path.hyPix)
				# Process climate
					∑Pet_Climate, ∑Pr_Climate, ∑T_Climate, N_∑T_Climate, Temp = climate.CLIMATE(clim, option)

			# OBSERVED θ  ~~~~~
			if option.hyPix.θobs
				# Read observed θ
					obsTheta = reading.hyPix.TIME_SERIES(clim, option, param, path.hyPix)
				
				# Process observed θ
					obsTheta = thetaObs.ΘOBS(obsTheta, clim, discret, Z)
			end #  option.hyPix.θobs

			# # MEMORY
			# 	if iOpt_Count == 1
			# 		Nit_Reduced = param.hyPix.iOpt_End - param.hyPix.iOpt_Start + 1

			# 		Efficiency                 = fill(0.0::Float64, Nit_Reduced)
			# 		Global_WaterBalance        = fill(0.0::Float64, Nit_Reduced)
			# 		Global_WaterBalance_NormPr = fill(0.0::Float64, Nit_Reduced)
			# 		RmseBest                   = fill(0.0::Float64, Nit_Reduced)
			# 		SwcRoots                   = fill(0.0::Float64, Nit_Reduced)
			# 		WofBest                    = fill(0.0::Float64, Nit_Reduced)
			# 		ΔRunTimeHypix              = fill(0.0::Float64, Nit_Reduced)
			# 		ΔT_Average                 = fill(0.0::Float64, Nit_Reduced)
			# 		∑ΔQ_Bot                    = fill(0.0::Float64, Nit_Reduced)
			# 		∑∑ΔSink                    = fill(0.0::Float64, Nit_Reduced)
			# 	end
				
				∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pr, ∑T, CropCoeficientᵀ, iNonConverge_iOpt, Laiᵀ, Q, Residual, ΔEvaporation, ΔHpond, ΔLnΨmax, ΔPet, ΔPr, ΔSink, ΔT, θ, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest = memory.MEMORY(clim, iOpt_Count, N_∑T_Climate, NiZ, obsTheta, param)

			# INITIALIZING THE STRUCTURE
				# Initializing hydroHorizon structure
					hydroHorizon = hydroStruct.HYDROSTRUCT(option.hyPix, N_Layer)

				# Initializing hydraulic param into structure 
					hydro = hydroStruct.HYDROSTRUCT(option.hyPix, NiZ)

				# Initialiozing vegetation parameters into veg structure
					veg = vegStruct.VEGSTRUCT()

				# Optimisation
					hydroHorizon_best = hydroStruct.HYDROSTRUCT(option.hyPix, N_Layer)
					hydro_best        = hydroStruct.HYDROSTRUCT(option.hyPix, NiZ)
					veg_best          = vegStruct.VEGSTRUCT()
				
			# ========================================================
	 
			# OBTAINING HYDRAULIC AND VEGETATION PARAMETERS
			if option.hyPix.Optimisation
				hydro, hydroHorizon, optim, veg = reading.hyPix.HYPIX_PARAM(Layer, hydro, hydroHorizon, iOpt, NiZ, option, param, path.hyPix.HyPixParamOpt, veg)
			else
			# Reading veg parameters
				veg, ~ = reading.READ_STRUCT(veg, path.hyPix.HyPix_VegParam)
	
			# Hydraulic parameters
				hydroHorizon, ~ = reading.READ_STRUCT(hydroHorizon, path.hyPix.HyPix_HydroParam)

				hydro           = horizonLayer.HYDROHORIZON_2_HYDRO(hydroHorizon, Layer, NiZ, option)

			# options of optim		
				Flag_Opt = false
				NparamOpt = 0

				optim = ( NparamOpt=NparamOpt, Flag_Opt=Flag_Opt)		
			end # option.hyPix.Optimisation

			# SINK TERM 
            Laiᵀ_η            = reading.hyPix.LOOKUPTABLE_LAI(clim, option, path.hyPix, veg)
            CropCoeficientᵀ_η = reading.hyPix.LOOKUPTABLE_CROPCOEFICIENT(clim, option, path.hyPix, veg)

			if optim.Flag_Opt
				hydro, hydro_best, hydroHorizon, hydroHorizon_best, veg, veg_best, WofBest = hypixOpt.HYPIXOPTIMISATION_START(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, hydro, hydro_best, hydroHorizon, hydroHorizon_best, iOpt_Count, Laiᵀ, Laiᵀ_η, Layer, N_∑T_Climate, N_Layer, NiZ, obsTheta, optim, option, param, Q, Residual, veg, veg_best, WofBest, Z, ΔEvaporation, ΔHpond, ΔPet, ΔPr, ΔSink, ΔT, ΔLnΨmax, θ, θini, θSim, Ψ, Ψini, Ψ_Max, Ψ_Min, Ψbest)
			end
			
			# if Flag_Opt then it will rerun with the optimal parameters
			println("\n =============== Start running HyPix =========================== \n")
			∑Pet, ∑Pr, ∑T, ∑T_Climate, clim, discret, iNonConverge, IterCount, N_iRoot, Nit, NiZ, Q, veg, ΔEvaporation, ΔHpond, ΔRootDensity, ΔT, θ, Ψ = hypixModel.HYPIX(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, hydro, Laiᵀ, Laiᵀ_η, N_∑T_Climate, NiZ, option, param, Q, Residual, veg, Z, ΔEvaporation, ΔHpond, ΔPet, ΔPr, ΔSink, ΔT, ΔLnΨmax, θ, θini, Ψ, Ψini, Ψ_Max, Ψ_Min, Ψbest)

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
					for iT=1:Nit
						∑ΔQ_Bot[iOpt_Count] = ∑ΔQ_Bot[iOpt_Count] + ΔT[iT] * Q[iT, NiZ+1]
					end
           
				println("		=== ===START: summary  $iOpt steps ...")

					println("			∑Pr 			= ", ceil(Int, ∑Pr_Clim), "  [mm]")
					println("			∑Pr_Soil 		= ", ceil(Int, ∑Pr[Nit] - ∑Pr[i∑T_CalibrStart]),  "  [mm]")
					# println("			∑Pr_Intercepted/∑Pr 	= ", ceil(Int, 100. * (∑Pr_Clim - (∑Pr[Nit]-∑Pr[i∑T_CalibrStart])) / (∑Pr_Clim + eps(10.0))),  "  [%]")
					println("			∑Pet_Net 		= ", ceil(Int, ∑Pet_Net), "  [mm]")
					println("			∑Pet 			= ", ceil(Int, ∑Pet[Nit]- ∑Pet[i∑T_CalibrStart]), "  [mm]")
					println("			∑ΔSink/∑Pet_Net 	= ", ceil(Int, 100.0 * ∑ΔSink[Nit] /(∑Pet_Net + eps(10.0))), "  [%] \n")
					
					println("			∑SoilWaterContentRootEnd = ", round(SwcRoots[iOpt_Count], digits=3), "  [mm]")
					println("			∑ΔSink 			= ", -ceil(Int, ∑∑ΔSink[iOpt_Count]), "  [mm]")
					println("			∑Infilt_Bot 		= ", -round(∑ΔQ_Bot[iOpt_Count],  digits=5), "  [mm]")
					println("			ΔHpond at end 		= ", ceil(Int, ΔHpond[Nit]), "  [mm] \n")

					println("			iNonConverge 			= ", iNonConverge_iOpt[iOpt_Count], "  [count]")
					println("			Global_WaterBalance_NormPr 	= ", round(Global_WaterBalance_NormPr[iOpt_Count], digits=8), "  [%]")
					println("			Global_WaterBalance 		= ", 	round(Global_WaterBalance[iOpt_Count], digits=8), "  [mm]")
					println("			Average ΔT 			= ",  ΔT_Average[iOpt_Count] , "  [seconds]")
					println("			ΔTmin 				= ",   round(minimum(ΔT[i∑T_CalibrStart:Nit]), digits=0) , "  [seconds]")
					println("			ΔTmax 				= ",  round(maximum(ΔT[i∑T_CalibrStart:Nit]), digits=0) , "  [seconds]")
					println("			ΔT_HyPix 			= ", ceil(Int, ΔRunTimeHypix[iOpt_Count]) , "  [seconds]")			
					println("			Efficiency 			= ", Efficiency[iOpt_Count], "  [iTer day-1]")
					println("			Number_of_cells 	    = ", NiZ, "  [-], \n")

					∑T_Reduced, ∑WaterBalanceη_Reduced, Date_Reduced, Nit_Reduced, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔPrGross_Reduced, ΔQ_Reduced, ΔSink_Reduced, ΔT_Reduced, θ_Reduced, θobs_Reduced, Ψ_Reduced = Δtchange.CHANGE_OUTPUT_ΔT(∑Pet[1:Nit], ∑Pr[1:Nit], ∑T[1:Nit], ∑WaterBalance_η[1:Nit], ∑ΔSink[1:Nit], obsTheta, clim, Nit, NiZ, param, Q[1:Nit,1:NiZ+1], ΔEvaporation[1:Nit], ΔHpond[1:Nit], ΔT[1:Nit], θ[1:Nit,1:NiZ], Ψ[1:Nit,1:NiZ], ∑T_Climate)

			# Computing average simulated θ to comapre it with average observed θ
			if option.hyPix.θobs_Average && option.hyPix.θobs	
				θsim_Aver = θaver.θAVER(discret; Z=Z, θ_Reduced=θ_Reduced, NiZ=NiZ, Nit_Reduced=Nit_Reduced, Zaver=min(400.0, Z[NiZ]))

				RmseBest[iOpt_Count] = stats.NSE(θobs_Reduced[1:Nit_Reduced], θsim_Aver[1:Nit_Reduced])

			elseif !(option.hyPix.θobs_Average) && option.hyPix.θobs
				RmseBest[iOpt_Count] = ofHypix.θof.RMSE_θ(∑T, obsTheta, Nit, NiZ, θ, θSim)
				θsim_Aver = []
				
			end

			if option.hyPix.θobs
				println("			Nse 			= ", round(RmseBest[iOpt_Count], digits=5), "  [mm3 mm-3]")
			end
			println("		=== === END: summary \n")


			if option.hyPix.Table
			println("		=== === START: Table === ===")
	
				# Writing values of hydraulic parameters
				table.hyPix.HYDRO(hydroHorizon, iOpt, N_Layer, path.hyPix)

				# Writing values of veg parameters
				table.hyPix.VEG(veg, iOpt, path.hyPix)

				table.hyPix.PERFORMANCE(∑∑ΔSink, ∑ΔQ_Bot, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iOpt, iOpt, RmseBest, SwcRoots, WofBest, ΔRunTimeHypix, ΔT_Average, Soilname, path.hyPix)		

				if option.hyPix.Table_Discretization
					table.hyPix.DISCRETIZATION(discret, NiZ, Z[1:NiZ], path.hyPix)
				end
				if  option.hyPix.Table_TimeSeries
					table.hyPix.TIME_SERIES(∑T[1:Nit], ΔT[1:Nit], ∑Pr[1:Nit], ΔPr[1:Nit], ΔHpond[1:Nit], ΔT[1:Nit].*Q[1:Nit,1], ∑WaterBalance_η[1:Nit], iOpt, path.hyPix)
				end
				if option.hyPix.Table_TimeSeriesDaily
					table.hyPix.TIME_SERIES_DAILY(∑T_Reduced[1:Nit_Reduced], ∑WaterBalanceη_Reduced[1:Nit_Reduced], Date_Reduced[1:Nit_Reduced], iOpt, Nit_Reduced, ΔEvaporation_Reduced[1:Nit_Reduced], ΔQ_Reduced[1:Nit_Reduced, NiZ+1], ΔPet_Reduced[1:Nit_Reduced], ΔPond_Reduced[1:Nit_Reduced], ΔPr_Reduced[1:Nit_Reduced], ΔSink_Reduced[1:Nit_Reduced], path.hyPix)
				end
				if option.hyPix.Table_θ
					table.hyPix.θ(Date_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced,1:NiZ], discret.Znode[1:NiZ], iOpt, path.hyPix)
				end
				if option.hyPix.Table_Ψ
					table.hyPix.Ψ(Date_Reduced[1:Nit_Reduced], Ψ_Reduced[1:Nit_Reduced,1:NiZ], discret.Znode[1:NiZ], iOpt, path.hyPix)
				end
				if option.hyPix.Table_Q
					table.hyPix.Q(Date_Reduced[1:Nit_Reduced], ΔQ_Reduced[1:Nit_Reduced,1:NiZ+1], Z[NiZ], discret.Znode[1:NiZ], iOpt, path.hyPix)
				end
				if option.hyPix.Tabule_θΨ
					table.hyPix.θΨ(hydroHorizon, iOpt, N_Layer, option.hyPix, param, path.hyPix)
					table.hyPix.KΨ(hydroHorizon, iOpt, N_Layer, option.hyPix, param, path.hyPix)
				end
				if option.hyPix.Table_Climate
					table.hyPix.DAILY_CLIMATE(∑T_Climate, clim, iOpt, path.hyPix)
				end
				if option.hyPix.θobs_Average && option.hyPix.θobs
					table.hyPix.θAVERAGE(Date_Reduced[1:Nit_Reduced], iOpt, θobs_Reduced[1:Nit_Reduced], θsim_Aver[1:Nit_Reduced], path.hyPix)
				end
			println("		=== === END: Table === === \n")
			end  # if option.hyPix.Table

		
			if option.other.Ploting
			println("		=== === START: Plotting === ===")

				# if option.hyPix.Plot_Other
				# plotOther.PLOT_θΨ_Δθ(hydro, path.hyPix, param, option)
				# 	# plotOther.ΨMINΨMAX(hydro, path.hyPix)
				# 	plotOther.WOF_STEPS(path.hyPix)
				# 	# plotOther.SE_Ψ_CONSTRAINED(hydro, path.hyPix)
				# 	# plotOther.PLOT_σ_2_θr(hydro, path.hyPix)
				# 	# plotOther.PLOT_θΨ_Δθ(hydro, path.hyPix)
				# end # option.hyPix.Plot_Other

			

				if option.hyPix.Plot_Hypix
					plotHypix.makkie.TIMESERIES(Date_Reduced, ∑T_Reduced, obsTheta, discret, iOpt, Nit_Reduced, NiZ, option, param, ΔEvaporation_Reduced, ΔQ_Reduced, ΔPrGross_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔSink_Reduced, θ_Reduced, θobs_Reduced, clim, i∑T_CalibrStart_Day, θsim_Aver, path.hyPix)

					if option.hyPix.Plot_θprofile
						plotHypix.makkie.θPROFILE(∑T_Reduced, discret, iSim, NiZ, obsTheta, option, param, path.hyPix, Soilname, θ_Reduced)
					end  # if: option.hyPix.Plot_
				end

				if option.hyPix.Plot_θΨK
					plotHypix.θΨK(hydroHorizon, N_Layer, iOpt, path.hyPix)
				end
				if option.hyPix.Plot_Vegetation && option.hyPix.RootWaterUptake
					plotHypix.VEG_FUNCTIONS(discret, iOpt, N_iRoot, veg, Z, ΔRootDensity, path.hyPix)
				end
				if option.hyPix.Plot_Interception
					plotHypix.plots.RAINFALL_INTERCEPTION(clim, i∑T_CalibrStart_Day, iOpt, path.hyPix)
				end
				if  option.hyPix.Plot_Sorptivity
					plotHypix.plots.PLOT_SORPTIVITY(hydro, iOpt, option, path.hyPix)
				end
			println("		=== === END: Plotting === === \n")
			end # if option.hyPix.Plotting
	
			println("	=== === === END  : Looping with time ")
			println("	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n")
		end # for loop: iOpt
	end #iSim = 1:length(Soilname)

	println("\n ==== END RUNNING Hypix_1D =====")

	return
	end # function HYPIX_START
	
end  # module hydro
# ............................................................