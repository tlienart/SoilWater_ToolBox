module psd
	import ..option, ..param, ..wrc, ..kunsat, ..cst, ..path, ..stats, ..psdFunc
	using Statistics, BlackBoxOptim, JuliaDB

	# ======================================================================================
	#          PSD_START
	# ======================================================================================
	function PSD_MAIN(N_SoilSelect, Ψ_θΨ, θ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, Rpart, ∑Psd, N_Psd, hydro)

		N_Psd_Max = maximum(N_Psd[1:N_SoilSelect])

		println("N_Psd_Max = $N_Psd_Max")
		println("option.psd.SubclayOpt = $(option.psd.SubclayOpt), option.psd.∑Psd_2_ξ2 = $(option.psd.∑Psd_2_ξ2), option.psd.∑Psd_2_ξ1 = $(option.psd.∑Psd_2_ξ1), \n")

		# plot.PLOTTING_MIXING_PARAM()
		# plot.PLOTTING_THETA_R_MODEL()
		# plot.PLOTTING_R_2_PSI_MODEL()
		# plot.PLOTTING_ξ2_MODEL()

		# # Number of soils wanted to run
		# 	N_SoilSelect = min(N_SoilSelect, param.psd.N_Soil_Select)

		# RESERVING MEMORY
			N_Psd 		= zeros(Int32, N_SoilSelect)

			Nse_Psd 	= zeros(Float64, N_SoilSelect)
			Nse_Psd_θh	= zeros(Float64, N_SoilSelect)
			Of_Psd		= zeros(Float64, N_SoilSelect)
			Subclay 	= zeros(Float64, N_SoilSelect)
			θr_Psd 		= zeros(Float64, N_SoilSelect)
			θr_Psd_Kg 	= zeros(Float64, N_SoilSelect)
			ξ1 			= zeros(Float64, N_SoilSelect)
			ξ2 			= zeros(Float64, N_SoilSelect)
		
			∑Psd 		= zeros(Float64, (N_SoilSelect, N_Psd_Max))
			Psd 		= zeros(Float64, (N_SoilSelect, N_Psd_Max))
			Rpart 		= zeros(Float64, (N_SoilSelect, N_Psd_Max))
			θ_Rpart 	= zeros(Float64, (N_SoilSelect, N_Psd_Max))
			Ψ_Rpart 	= zeros(Float64, (N_SoilSelect, N_Psd_Max))

			Rpart_Whole = zeros(Float64, N_Psd_Max)
			ξ 			= zeros(Float64, N_Psd_Max) # Problem need to solve

		# INITIAL VALUES WHICH CAN CHANGE (depending on options)
			∑Psd_2_ξ2_β1 = param.psd.∑Psd_2_ξ2_β1
			∑Psd_2_ξ2_β2 = param.psd.∑Psd_2_ξ2_β2

		# Computing the radius of the particles from the diameters
			# Rpart_Whole = param.psd.Dpart ./ 2.0

		# Reading PSD data
			# ∑Psd = reading.PSD(path.Psd, Flag_Good) # now need to read with the general reading function!!!!!!!!! 
			# Diameter_Psd, ∑Psd, N_Psd = read.PSD(Id_True, N_SoilSelect)

		# =================== FOR EVERY SOIL SAMPLES ===================
		for iSoil=1:N_SoilSelect
			# Compute maximum N particle size for each iSoil

			println(∑Psd[iSoil,1:N_Psd[iSoil]])

			# Correction for N_Psd such that to determine the real maximum Rpart size
			# N_Psd_New = 1
			# for iPsd = 1:N_Psd[iSoil]
			# 	# println("$iSoil, $iPsd, $(∑Psd[iSoil, iPsd])")
			# 	# if ∑Psd[iSoil, iPsd] >= 0.99999
			# 	# 	N_Psd_New = iPsd
			# 	# 	break
			# 	# end
			# end
			# N_Psd[iSoil] = N_Psd_New
			# println("$iSoil $(N_Psd[iSoil])")
			
			# Compute PSD from ∑PSD
			Psd[iSoil,1:N_Psd[iSoil]] = psdFunc.∑PSD_2_PSD(∑Psd[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil])

		end # looping over soils


		# =================== COMPUTTING  θr model from Psd ===================
			if option.psd.Psd_2_θr == "Opt" 
				θr_Psd = psdThetar.OPTIMIZE_PSD_2_θr(N_SoilSelect, θr[1:N_SoilSelect], ∑Psd[1:N_SoilSelect,:])

			elseif  option.psd.Psd_2_θr == "Cst"
				@simd for iSoil=1:N_SoilSelect
					θr_Psd[iSoil] = param.psd.θr_Cst
				end
				
			elseif  option.psd.Psd_2_θr == "Param"
				println("Optimize θr = Psd_2_θr_α1 = $(param.psd.Psd_2_θr_α1) ; param.psd.Psd_2_θr_α2 = $(param.psd.Psd_2_θr_α2)")
				@simd for iSoil=1:N_SoilSelect
					θr_Psd[iSoil] = psdFunc.PSD_2_θr(∑Psd[iSoil,param.psd.Psd_2_θr_Size])
				end
			else
				error("option.psd.Psd_2_θr = $option.psd.Psd_2_θr  not allowed option.psd.Psd_2_θr must be either (1)'Opt' or (2) 'Cst' or (3) 'Param' ")
			end # if option.psd.Psd_2_θr

			Nse_θr = 1.0 - stats.NASH_SUTCLIFE_MINIMIZE(θr[1:N_SoilSelect], θr_Psd[1:N_SoilSelect])

			println("NASH_SUTCLIFFE θr_Psd = $Nse_θr \n")

			if option.psd.Plot_θr
				plot.PLOTTING_θr(θr[1:N_SoilSelect], θr_Psd[1:N_SoilSelect], ∑Psd[1:N_SoilSelect,param.psd.Psd_2_θr_Size])
			end # option.psd.Plot_θr

		# ........................................................................


		if option.psd.OptimizePsd == "Single" # <> <> <> <> <> <>
			∑Of_Psd = 0.0

			@simd for iSoil=1:N_SoilSelect
				ξ1[iSoil], ξ2[iSoil], Ψ_Rpart[iSoil,1:N_Psd[iSoil]], θ_Rpart[iSoil,1:N_Psd[iSoil]], Of_Psd[iSoil], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Subclay[iSoil] = OPTIMIZATION_SINGLE_SOIL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil], θr_Psd[iSoil])
			end # Loop single

		elseif option.psd.OptimizePsd == "All" # <> <> <> <> <> <>
			∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd = OPTIMIZATION_ALL_SOIL(N_SoilSelect, Psd[1:N_SoilSelect,:], ∑Psd[1:N_SoilSelect,:], Rpart[1:N_SoilSelect,:], N_Psd[1:N_SoilSelect], θsMac[1:N_SoilSelect], θr[1:N_SoilSelect], ΨkgMat[1:N_SoilSelect], σMat[1:N_SoilSelect], θsMat[1:N_SoilSelect], ΨkgMac[1:N_SoilSelect], σMac[1:N_SoilSelect], θr_Psd[1:N_SoilSelect])
		
		elseif option.psd.OptimizePsd == "Run" # <> <> <> <> <> <>
			∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd = PARAMETERS_ALL_SOIL(N_SoilSelect, Psd[1:N_SoilSelect,:], ∑Psd[1:N_SoilSelect,:], Rpart[1:N_SoilSelect,:], N_Psd[1:N_SoilSelect], θsMac[1:N_SoilSelect], θr[1:N_SoilSelect], ΨkgMat[1:N_SoilSelect], σMat[1:N_SoilSelect], θsMat[1:N_SoilSelect], ΨkgMac[1:N_SoilSelect], σMac[1:N_SoilSelect], θr_Psd[1:N_SoilSelect])

		else
			error("  $(option.psd.OptimizePsd) not found ")
		end # option.psd.OptimizePsd

		# STATISTICS NSE
			println("\n PSD MODEL:")
			Nse_Psd, Nse_Psd_Mean = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd[1:N_SoilSelect], Ψ_Rpart[1:N_SoilSelect,:], θ_Rpart[1:N_SoilSelect,:], θsMac[1:N_SoilSelect], θr[1:N_SoilSelect], ΨkgMat[1:N_SoilSelect], σMat[1:N_SoilSelect], θsMat[1:N_SoilSelect], ΨkgMac[1:N_SoilSelect], σMac[1:N_SoilSelect])

			# Nse_Psd_Mean = Statistics.mean(Nse_Psd[1:N_SoilSelect])


		# =================================================================================================================
		#       OPTIMIZATION HYDRAULIC PARAMETERS DERIVED FROM psd
		# =================================================================================================================
		if option.psd.HydroParam
			println("START Optimizing the hydraulic parameters derived from PSD")
			println("...")

			θsMat_Psd, θr_Psd_Kg, σMat_Psd, ΨkgMat_Psd, ~, θsMac_Psd, σMac_Psd, ΨkgMac_Psd, ~, Nse_θh_Uni_Psd, Nse_θh_Bim_Psd = optimiseCharacUnsat.OPTIMISE_CHARAC_UNSAT(Ψ_Rpart[1:N_SoilSelect,:], θ_Rpart[1:N_SoilSelect,:], N_Psd[1:N_SoilSelect], N_SoilSelect, θsMac[1:N_SoilSelect], θr_Psd=θr_Psd[1:N_SoilSelect]; Option_Data_Kθ = false )

			Nse_Psd_σ = 1. - stats.NASH_SUTCLIFE_MINIMIZE(σMat, σMat_Psd; Power=2.)
			Nse_Psd_Ψkg = 1. - stats.NASH_SUTCLIFE_MINIMIZE(ΨkgMat, ΨkgMat_Psd; Power=2.)
			ΔθsMacMat = θsMac[1:N_SoilSelect] .- θsMat[1:N_SoilSelect] 
			ΔθsMacMat_Psd = θsMac[1:N_SoilSelect] .- θsMat_Psd[1:N_SoilSelect] 
			Nse_Psd_θs = 1. - stats.NASH_SUTCLIFE_MINIMIZE(ΔθsMacMat, ΔθsMacMat_Psd ; Power=2.)

			println("Nse_Psd_σ = $Nse_Psd_σ ; Nse_Psd_Ψkg = $Nse_Psd_Ψkg ; Nse_Psd_θs =$Nse_Psd_θs")

			table.SIMULATION(Nse_Psd_σ, Nse_Psd_Ψkg, Nse_Psd_θs, Nse_Psd_Mean, Nse_θΨ_Bim_Mean)

			if  option.psd.Plot_σ_Ψkg
				plot.PLOTTING_σ_Ψkg(σMat[1:N_SoilSelect], σMat_Psd[1:N_SoilSelect], ΨkgMat[1:N_SoilSelect], ΨkgMat_Psd[1:N_SoilSelect], ΔθsMacMat[1:N_SoilSelect], ΔθsMacMat_Psd[1:N_SoilSelect])
			end # option.psd.Plot_σ_Ψkg

			println("END Optimizing the hydraulic parameters derived from PSD \n") 
		end # option.Psd.HydroParam
	

		# =================================================================================================================
		#			TABLES
		# =================================================================================================================
			println("START WRITTING TABLE, \n")
			println(" ... ")

			if option.psd.OptimizePsd == "Single" 
				table.SINGLEOPT_T1_T2(θsMac[1:N_SoilSelect], θr[1:N_SoilSelect], θr_Psd[1:N_SoilSelect], σMat[1:N_SoilSelect], ΨkgMat[1:N_SoilSelect], θsMat[1:N_SoilSelect], σMac[1:N_SoilSelect], ΨkgMac[1:N_SoilSelect], ξ1[1:N_SoilSelect], ξ2[1:N_SoilSelect], Nse_Psd[1:N_SoilSelect], Subclay[1:N_SoilSelect])
			end

			if option.psd.HydroParam
				table.HYDRAULICPARAM_Psd(θsMat_Psd[1:N_SoilSelect], θr_Psd_Kg[1:N_SoilSelect], σMat_Psd[1:N_SoilSelect], ΨkgMat_Psd[1:N_SoilSelect], θsMac_Psd[1:N_SoilSelect], σMac_Psd[1:N_SoilSelect], ΨkgMac_Psd[1:N_SoilSelect], Nse_θh_Uni_Psd[1:N_SoilSelect],Nse_θh_Bim_Psd[1:N_SoilSelect], ∑Psd[1:N_SoilSelect,1:param.psd.N_Psd])
			end

			table.INTERGRANULARMIXING(θsMac[1:N_SoilSelect], θr[1:N_SoilSelect], θr_Psd[1:N_SoilSelect], σMat[1:N_SoilSelect], ΨkgMat[1:N_SoilSelect], θsMat[1:N_SoilSelect], σMac[1:N_SoilSelect], ΨkgMac[1:N_SoilSelect], ξ1[1:N_SoilSelect], ξ2[1:N_SoilSelect], Nse_Psd[1:N_SoilSelect], Subclay[1:N_SoilSelect], ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, 0.0, 0.0, 0.0)

			println("END WRITTING TABLE, \n")


		# =================================================================================================================
		#			PLOTTING
		# =================================================================================================================
			
		plot.PLOTTING_MIXING_PARAM(;ξ1=6.56, ξ2=0.18)
		
			if option.psd.Plot_∑Psd_2_ξ2
				Data = JuliaDB.loadtable(path.Table * "\\IndividualSamples_T1_T2.csv")
				ξ2_Obs = JuliaDB.select(Data,:Tau2)
				plot.PLOTTING_∑Psd_2_ξ2(ξ2_Obs[1:N_SoilSelect], ξ2[1:N_SoilSelect], ∑Psd[1:N_SoilSelect, param.psd.∑Psd_2_ξ2_Size])

				#  Paper graph
				plot.PLOTTING_ξ2_MODEL(ξ2_Obs[1:N_SoilSelect], ξ2[1:N_SoilSelect], ∑Psd[1:N_SoilSelect, param.psd.∑Psd_2_ξ2_Size])

				# Statistics
				Nse_ξ2 = 1.0 - stats.NASH_SUTCLIFE_MINIMIZE(ξ2_Obs[1:N_SoilSelect], ξ2[1:N_SoilSelect])
				println("NSE_ξ2 = $Nse_ξ2")

			end # option.psd.Plot_∑Psd_2_ξ2

			if option.psd.Plot_OFsingle_OFmeas
				Data = JuliaDB.loadtable(path.Table * "\\IndividualSamples_T1_T2.csv")
				Nse_T1T2 = JuliaDB.select(Data,:NSE_Psd)

				Data = JuliaDB.loadtable(path.OptimizeKΨ)
				Nse_Obs = JuliaDB.select(Data,:Nse_THETAh_Bim)

				plot.PLOTTING_OF(Nse_T1T2[1:N_SoilSelect], Nse_Obs[1:N_SoilSelect], Nse_Psd[1:N_SoilSelect])
			end # option.psd.Plot_OFsingle_OFmeas

			if option.psd.Plotting
				@simd for iSoil = 1:N_SoilSelect
					println("Plotting soil $iSoil")

					#   INTERGRANULARMIXING
					@simd for iRpart = 1: N_Psd[iSoil]
						ξ[iRpart] = psdFunc.INTERGRANULARMIXING(Rpart[iSoil,iRpart], ξ1[iSoil], ξ2[iSoil])
					end

					if option.psd.OptimizeKΨ && option.psd.HydroParam
						plot.PLOTTING(iSoil, N_Psd[iSoil], Rpart[iSoil,1:N_Psd[iSoil]], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], θsMac[iSoil], θr[iSoil], σMat[iSoil], ΨkgMat[iSoil], KsMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMac[iSoil], ΨkgMac[iSoil], Ψ_θΨ[iSoil,1:N_θΨ[iSoil]], θ_θΨ[iSoil,1:N_θΨ[iSoil]], N_θΨ[iSoil], Ψ_Rpart[iSoil,1:N_Psd[iSoil]], θ_Rpart[iSoil,1:N_Psd[iSoil]], ξ[1:N_Psd[iSoil]]; N_Kθ=N_Kθ[iSoil], K_Kθ=K_Kθ[iSoil,1:N_Kθ[iSoil]], Ψ_Kθ=Ψ_Kθ[iSoil,1:N_Kθ[iSoil]],θsMat_Psd=θsMat_Psd[iSoil], θr_Psd=θr_Psd[iSoil], σMat_Psd=σMat_Psd[iSoil], ΨkgMat_Psd=ΨkgMat_Psd[iSoil], θsMac_Psd=θsMac_Psd[iSoil], σMac_Psd=σMac_Psd[iSoil], ΨkgMac_Psd=ΨkgMac_Psd[iSoil])

					elseif option.psd.OptimizeKΨ && !option.psd.HydroParam
						plot.PLOTTING(iSoil, N_Psd[iSoil], Rpart[iSoil,1:N_Psd[iSoil]], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], θsMac[iSoil], θr[iSoil], σMat[iSoil], ΨkgMat[iSoil], KsMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMac[iSoil], ΨkgMac[iSoil], Ψ_θΨ[iSoil,1:N_θΨ[iSoil]], θ_θΨ[iSoil,1:N_θΨ[iSoil]], N_θΨ[iSoil], Ψ_Rpart[iSoil,1:N_Psd[iSoil]], θ_Rpart[iSoil,1:N_Psd[iSoil]], ξ[1:N_Psd[iSoil]]; N_Kθ=N_Kθ[iSoil], K_Kθ=K_Kθ[iSoil,1:N_Kθ[iSoil]], Ψ_Kθ=Ψ_Kθ[iSoil,1:N_Kθ[iSoil]])

					elseif !option.psd.OptimizeKΨ && option.psd.HydroParam
						plot.PLOTTING(iSoil, N_Psd[iSoil], Rpart[iSoil,1:N_Psd[iSoil]],Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], θsMac[iSoil], θr[iSoil], σMat[iSoil], ΨkgMat[iSoil], KsMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMac[iSoil], ΨkgMac[iSoil], Ψ_θΨ[iSoil,1:N_θΨ[iSoil]], θ_θΨ[iSoil,1:N_θΨ[iSoil]], N_θΨ[iSoil], Ψ_Rpart[iSoil,1:N_Psd[iSoil]], θ_Rpart[iSoil,1:N_Psd[iSoil]], ξ[1:N_Psd[iSoil]];θsMat_Psd=θsMat_Psd[iSoil], θr_Psd=θr_Psd[iSoil], σMat_Psd=σMat_Psd[iSoil], ΨkgMat_Psd=ΨkgMat_Psd[iSoil], θsMac_Psd=θsMac_Psd[iSoil], σMac_Psd=σMac_Psd[iSoil], ΨkgMac_Psd=ΨkgMac_Psd[iSoil])
					
					elseif !option.psd.OptimizeKΨ && !option.psd.HydroParam
						plot.PLOTTING(iSoil, N_Psd[iSoil], Rpart[iSoil,1:N_Psd[iSoil]], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], θsMac[iSoil], θr[iSoil], σMat[iSoil], ΨkgMat[iSoil], KsMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMac[iSoil], ΨkgMac[iSoil], Ψ_θΨ[iSoil,1:N_θΨ[iSoil]], θ_θΨ[iSoil,1:N_θΨ[iSoil]], N_θΨ[iSoil], Ψ_Rpart[iSoil,1:N_Psd[iSoil]], θ_Rpart[iSoil,1:N_Psd[iSoil]], ξ[1:N_Psd[iSoil]])
					end # elseif

				end # looping over soils
			end #  if option.psd.Plotting
		return
	end # function PSD_START



	# =================================================================================================================
	#       OPTIMIZATION SINGLE SOIL
	# =================================================================================================================
	function OPTIMIZATION_SINGLE_SOIL(iSoil, Psd, ∑Psd, Rpart, N_Psd, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)
		if option.psd.SubclayOpt  && !(option.psd.∑Psd_2_ξ2) && !(option.psd.∑Psd_2_ξ1) # <><><><><><><>	

			Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(P[1], P[2], P[3], Rpart[1:N_Psd], N_Psd, Psd[1:N_Psd], ∑Psd[1:N_Psd], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd) ; SearchRange =[(param.psd.ξ1_Min, param.psd.ξ1_Max), (param.psd.ξ2_Min, param.psd.ξ2_Max), (param.psd.Wsubclay_Min, param.psd.Wsubclay_Max)], NumDimensions=3,  TraceMode=:silent)

			 # Optimal values
			ξ1 = BlackBoxOptim.best_candidate(Optimization)[1]
			ξ2 = BlackBoxOptim.best_candidate(Optimization)[2]
			Subclay = BlackBoxOptim.best_candidate(Optimization)[3]
			Of_Psd = BlackBoxOptim.best_fitness(Optimization)

		elseif !(option.psd.SubclayOpt) && !(option.psd.∑Psd_2_ξ2) && !(option.psd.∑Psd_2_ξ1)# <><><><><><><># <><><><><><><>
			Subclay = param.psd.Subclay
	
			Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(P[1], P[2], Subclay , Rpart[1:N_Psd], N_Psd, Psd[1:N_Psd], ∑Psd[1:N_Psd], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd) ; SearchRange =[(param.psd.ξ1_Min, param.psd.ξ1_Max), (param.psd.ξ2_Min, param.psd.ξ2_Max)], NumDimensions=2,  TraceMode=:silent)

			 # Optimal values
			ξ1 = BlackBoxOptim.best_candidate(Optimization)[1]
			ξ2 = BlackBoxOptim.best_candidate(Optimization)[2]
			Of_Psd = BlackBoxOptim.best_fitness(Optimization)

		elseif !(option.psd.SubclayOpt) && !(option.psd.∑Psd_2_ξ2)  && (option.psd.∑Psd_2_ξ1)
			Subclay = param.psd.Subclay

			ξ1 = param.psd.P_ξ1
	
			Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(ξ1 , P[1], Subclay , Rpart[1:N_Psd], N_Psd, Psd[1:N_Psd], ∑Psd[1:N_Psd], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd) ; SearchRange =[(param.psd.ξ2_Min, param.psd.ξ2_Max)], NumDimensions=1,  TraceMode=:silent)
			 # Optimal values
			ξ2 = BlackBoxOptim.best_candidate(Optimization)[1]
			Of_Psd = BlackBoxOptim.best_fitness(Optimization)


		elseif option.psd.∑Psd_2_ξ2  && !(option.psd.∑Psd_2_ξ1)
			Subclay = param.psd.Subclay
			ξ2 = psdFunc.∑PSD_2_ξ2(∑Psd[param.psd.∑Psd_2_ξ2_Size])

			Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(P[1], ξ2 , Subclay , Rpart[1:N_Psd], N_Psd, Psd[1:N_Psd], ∑Psd[1:N_Psd], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd) ; SearchRange =[(param.psd.ξ1_Min, param.psd.ξ1_Max)], NumDimensions=1,  TraceMode=:silent)

			 # Optimal values
			ξ1 = BlackBoxOptim.best_candidate(Optimization)[1]
			Of_Psd = BlackBoxOptim.best_fitness(Optimization)


		elseif option.psd.∑Psd_2_ξ2  && option.psd.∑Psd_2_ξ1
			Subclay = param.psd.Subclay # Average value

			ξ2 =  psdFunc.∑PSD_2_ξ2(∑Psd[param.psd.∑Psd_2_ξ2_Size])

			ξ1 = param.psd.P_ξ1

			Of_Psd = OF_SINGLE_SOIL(ξ1, ξ2, param.psd.Subclay, Rpart[1:N_Psd], N_Psd, Psd[1:N_Psd], ∑Psd[1:N_Psd], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)
		else
			error(" NO SUITABLE OPTIONS FOUND ")
		end # if option.psd.

		# Recording the optimal parameters
		θ_Rpart, Ψ_Rpart = psdFunc.PSD_MODEL(Rpart, Psd, ∑Psd, N_Psd, θsMac, θr_Psd, Subclay, ξ1, ξ2)

		# Correction of Subclay  and deriving the corrected Psd and  ∑Psd
		Psd, ∑Psd = psdFunc.SUBCLAY_CORRECTION(∑Psd, Subclay, N_Psd)

		return ξ1, ξ2, Ψ_Rpart, θ_Rpart, Of_Psd, Psd, ∑Psd, Subclay
	end # OPTIMIZATION_SINGLE_SOIL



	function OF_SINGLE_SOIL(ξ1, ξ2, Subclay, Rpart, N_Psd, Psd, ∑Psd, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd) 
		θ_Rpart = zeros(Float64, N_Psd)
		Ψ_Rpart = zeros(Float64, N_Psd)
		θΨ = zeros(Float64, N_Psd)

		θ_Rpart[1:N_Psd], Ψ_Rpart[1:N_Psd] = psdFunc.PSD_MODEL(Rpart, Psd, ∑Psd, N_Psd, θsMac, θr_Psd, Subclay, ξ1, ξ2)

		# For every class
		Of = 0.0
		for iRpart = 1:N_Psd
			# Observed data
			θΨ[iRpart] = wrc.kg.Ψ_2_θdual(Ψ_Rpart[iRpart], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac)
			# Summing the error
			Of += (θΨ[iRpart] - θ_Rpart[iRpart])^2.0
		end
		Of = (Of / N_Psd) ^ 0.5
		return Of
	end # OF_SINGLE_SOIL ===============



	#= =================================================================================================================
	       Optimisation ALL SOILS
	================================================================================================================= =#
	 function OPTIMIZATION_ALL_SOIL(N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)

		function OF_ALL_SOIL(ξ1, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)
			Of_AllSoil = 0.
			@simd for iSoil = 1:N_SoilSelect
				ξ2 =  psdFunc.∑PSD_2_ξ2(∑Psd[iSoil, param.psd.∑Psd_2_ξ2_Size]; ∑Psd_2_ξ2_β1=∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=∑Psd_2_ξ2_β2) #####################for new table model 4######
				# ξ2 = ∑Psd_2_ξ2_β1   ############################for new table model 1######

				Of_AllSoil += OF_SINGLE_SOIL(ξ1, ξ2, Subclay, Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil], θr_Psd[iSoil])
			end # for loop
			return Of_AllSoil
		end # function OF_ALL_SOILS


		if option.psd.SubclayOpt 
			Optimization = BlackBoxOptim.bboptimize(P -> OF_ALL_SOIL(P[1], P[2], P[3], P[4], N_SoilSelect, Psd[1:N_SoilSelect,:], ∑Psd[1:N_SoilSelect,:], Rpart[1:N_SoilSelect,:], N_Psd[1:N_SoilSelect], θsMac[1:N_SoilSelect], θr[1:N_SoilSelect], ΨkgMat[1:N_SoilSelect], σMat[1:N_SoilSelect], θsMat[1:N_SoilSelect], ΨkgMac[1:N_SoilSelect], σMac[1:N_SoilSelect], θr_Psd[1:N_SoilSelect]);  SearchRange =[(param.psd.ξ1_Min, param.psd.ξ1_Max), (param.psd.∑Psd_2_ξ2_β1_Min,param.psd.∑Psd_2_ξ2_β1_Max), (param.psd.∑Psd_2_ξ2_β2_Min,param.psd.∑Psd_2_ξ2_β2_Max), (param.psd.Wsubclay_Min, param.psd.Wsubclay_Max)], NumDimensions=4, TraceMode=:silent )

			ξ1_All =  BlackBoxOptim.best_candidate(Optimization)[1]
			∑Psd_2_ξ2_β1 = BlackBoxOptim.best_candidate(Optimization)[2]
			∑Psd_2_ξ2_β2 =  BlackBoxOptim.best_candidate(Optimization)[3]
			Wsubclay_All =  BlackBoxOptim.best_candidate(Optimization)[4]

		elseif !(option.psd.SubclayOpt)
			Wsubclay_All = param.psd.Subclay

			Optimization = BlackBoxOptim.bboptimize(P -> OF_ALL_SOIL(P[1], P[2], P[3], Wsubclay_All, N_SoilSelect, Psd[1:N_SoilSelect,:], ∑Psd[1:N_SoilSelect,:], Rpart[1:N_SoilSelect,:], N_Psd[1:N_SoilSelect], θsMac[1:N_SoilSelect], θr[1:N_SoilSelect], ΨkgMat[1:N_SoilSelect], σMat[1:N_SoilSelect], θsMat[1:N_SoilSelect], ΨkgMac[1:N_SoilSelect], σMac[1:N_SoilSelect], θr_Psd[1:N_SoilSelect]);  SearchRange =[(param.psd.ξ1_Min, param.psd.ξ1_Max), (param.psd.∑Psd_2_ξ2_β1_Min,param.psd.∑Psd_2_ξ2_β1_Max), (param.psd.∑Psd_2_ξ2_β2_Min,param.psd.∑Psd_2_ξ2_β2_Max)], NumDimensions=3, TraceMode=:silent )

			ξ1_All =  BlackBoxOptim.best_candidate(Optimization)[1]
			∑Psd_2_ξ2_β1 = BlackBoxOptim.best_candidate(Optimization)[2]
			∑Psd_2_ξ2_β2 =  BlackBoxOptim.best_candidate(Optimization)[3]
		end # if option.psd.SubclayOpt 

		Of_AllSoil = BlackBoxOptim.best_fitness(Optimization)

		println("OPTIMAL PSD HYDRAULIC PARAMETERS:")
		println("ξ1=$ξ1_All")
		println("∑Psd_2_ξ2_β1=$∑Psd_2_ξ2_β1")
		println("∑Psd_2_ξ2_β2=$∑Psd_2_ξ2_β2")
		# println("NSE_ξ2 = $Nse_ξ2")          hay que introducir la variable en la funcion  !!!!!!!!
		println("Subclay=$Wsubclay_All")

		# Optimal parameters
		local N_Psd = 100 # Just a guess
		ξ1 = zeros(Float64, N_SoilSelect)
		ξ2 = zeros(Float64, N_SoilSelect)
		Subclay = zeros(Float64, N_SoilSelect)
		Ψ_Rpart = zeros(Float64, N_SoilSelect, N_Psd)
		θ_Rpart = zeros(Float64, N_SoilSelect, N_Psd)
		Of_Psd = zeros(Float64, N_SoilSelect, N_Psd)

		@simd for iSoil = 1:N_SoilSelect
			Subclay[iSoil] = Wsubclay_All

			ξ1[iSoil] =  ξ1_All 

			ξ2[iSoil] =  psdFunc.∑PSD_2_ξ2(∑Psd[iSoil, param.psd.∑Psd_2_ξ2_Size]; ∑Psd_2_ξ2_β1=∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=∑Psd_2_ξ2_β2) #####################for new table model 4######
			# ξ2[iSoil] = ∑Psd_2_ξ2_β1  #####################for new table model 1######

			θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] = psdFunc.PSD_MODEL(Rpart[iSoil,1:N_Psd[iSoil]], Psd[iSoil,:], ∑Psd[iSoil,:], N_Psd[iSoil], θsMac[iSoil], θr_Psd[iSoil], Subclay[iSoil], ξ1[iSoil], ξ2[iSoil])

			Of_Psd[iSoil] = OF_SINGLE_SOIL(ξ1[iSoil], ξ2[iSoil], Subclay[iSoil], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil], θr_Psd[iSoil])
		end # For

		return ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd
	end # function OPTIMIZATION_ALL_SOIL

#= =================================================================================================================
	       Parameters ALL SOILS
	================================================================================================================= =#
	function PARAMETERS_ALL_SOIL(N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)

		function OF_ALL_SOIL(ξ1, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)
			Of_AllSoil = 0.0
			@simd for iSoil = 1:N_SoilSelect
				ξ2 =  psdFunc.∑PSD_2_ξ2(∑Psd[iSoil, param.psd.∑Psd_2_ξ2_Size]; ∑Psd_2_ξ2_β1=param.psd.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=param.psd.∑Psd_2_ξ2_β2) 				

				Of_AllSoil += OF_SINGLE_SOIL(ξ1, ξ2, Subclay, Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil], θr_Psd[iSoil])
			end # for loop
			return Of_AllSoil
		end # function OF_ALL_SOILS


		if option.psd.SubclayOpt 
			
			ξ1_All =  param.psd.ξ1
			∑Psd_2_ξ2_β1 = param.psd.∑Psd_2_ξ2_β1
			∑Psd_2_ξ2_β2 = param.psd.∑Psd_2_ξ2_β2
			Wsubclay_All =  param.psd.Subclay

		elseif !(option.psd.SubclayOpt)
			Wsubclay_All = param.psd.Subclay

			ξ1_All =  param.psd.ξ1
			∑Psd_2_ξ2_β1 = param.psd.∑Psd_2_ξ2_β1
			∑Psd_2_ξ2_β2 = param.psd.∑Psd_2_ξ2_β2
		end # if option.psd.SubclayOpt 

		println("PSD HYDRAULIC PARAMETERS from PARAM:")
		println("ξ1=$ξ1_All")
		println("∑Psd_2_ξ2_β1=$∑Psd_2_ξ2_β1")
		println("∑Psd_2_ξ2_β2=$∑Psd_2_ξ2_β2")
		println("Subclay=$Wsubclay_All")

		# Optimal parameters
		local N_Psd = 100 # Just a guess
		ξ1 = zeros(Float64, N_SoilSelect)
		ξ2 = zeros(Float64, N_SoilSelect)
		Subclay = zeros(Float64, N_SoilSelect)
		Ψ_Rpart = zeros(Float64, N_SoilSelect, N_Psd)
		θ_Rpart = zeros(Float64, N_SoilSelect, N_Psd)
		Of_Psd = zeros(Float64, N_SoilSelect, N_Psd)

		@simd for iSoil = 1:N_SoilSelect
			Subclay[iSoil] = Wsubclay_All

			ξ1[iSoil] =  ξ1_All 

			ξ2[iSoil] =  psdFunc.∑PSD_2_ξ2(∑Psd[iSoil, param.psd.∑Psd_2_ξ2_Size]; ∑Psd_2_ξ2_β1=param.psd.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=param.psd.∑Psd_2_ξ2_β2) 			

			θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] = psdFunc.PSD_MODEL(Rpart[iSoil,1:N_Psd[iSoil]], Psd[iSoil,:], ∑Psd[iSoil,:], N_Psd[iSoil], θsMac[iSoil], θr_Psd[iSoil], Subclay[iSoil], ξ1[iSoil], ξ2[iSoil])

			Of_Psd[iSoil] = OF_SINGLE_SOIL(ξ1[iSoil], ξ2[iSoil], Subclay[iSoil], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil], θr_Psd[iSoil])
		end # For

		return ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd
	end # function PARAMETERS_ALL_SOIL

end # module PSD