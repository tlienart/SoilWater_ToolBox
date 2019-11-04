module psd
	import ..option, ..param, ..wrc, ..kunsat, ..cst, ..path, ..stats, ..psdFunc, ..psdInitiate, ..psdThetar, ..table
	using Statistics, BlackBoxOptim, JuliaDB

	# ======================================================================================
	#          PSD_START
	# ======================================================================================
	function PSD_MAIN(N_SoilSelect, Ψ_θΨ, θ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, Rpart, ∑Psd, N_Psd, hydro)

		# INITIATING THE PSD DATA		
			N_Psd, N_Psd_Max, Psd = psdInitiate.PSD_INITIATE(N_Psd, N_SoilSelect, ∑Psd)
		
		# COMPUTING θr FROM PSD DATA
			Nse_θr, θr_Psd = psd.psdThetar.MAIN_PSD_2_θr(N_SoilSelect, ∑Psd, hydro)
			println("θr_Psd = $θr_Psd")
			println("NASH_SUTCLIFFE θr_Psd = $Nse_θr \n")
		
		# 



		if option.psd.OptimizePsd == "Single" # <> <> <> <> <> <>
			∑Of_Psd = 0.0

			@simd for iSoil=1:N_SoilSelect
				ξ1[iSoil], ξ2[iSoil], Ψ_Rpart[iSoil,1:N_Psd[iSoil]], θ_Rpart[iSoil,1:N_Psd[iSoil]], Of_Psd[iSoil], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Subclay[iSoil] = OPTIMIZATION_SINGLE_SOIL(iSoil, Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], hydro.θsMac[iSoil], hydro.θr[iSoil], hydro.ΨkgMat[iSoil], hydro.σMat[iSoil], hydro.θsMat[iSoil], hydro.ΨkgMac[iSoil], hydro.σMac[iSoil], θr_Psd[iSoil])
			end # Loop single

		elseif option.psd.OptimizePsd == "All" # <> <> <> <> <> <>
			∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd = OPTIMIZATION_ALL_SOIL(N_SoilSelect, Psd[1:N_SoilSelect,:], ∑Psd[1:N_SoilSelect,:], Rpart[1:N_SoilSelect,:], N_Psd[1:N_SoilSelect], hydro.θsMac[1:N_SoilSelect], hydro.θr[1:N_SoilSelect], hydro.ΨkgMat[1:N_SoilSelect], hydro.σMat[1:N_SoilSelect], hydro.θsMat[1:N_SoilSelect], hydro.ΨkgMac[1:N_SoilSelect], hydro.σMac[1:N_SoilSelect], θr_Psd[1:N_SoilSelect])
		
		elseif option.psd.OptimizePsd == "Run" # <> <> <> <> <> <>
			∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd = PARAMETERS_ALL_SOIL(N_SoilSelect, Psd[1:N_SoilSelect,:], ∑Psd[1:N_SoilSelect,:], Rpart[1:N_SoilSelect,:], N_Psd[1:N_SoilSelect], hydro.θsMac[1:N_SoilSelect], hydro.θr[1:N_SoilSelect], hydro.ΨkgMat[1:N_SoilSelect], hydro.σMat[1:N_SoilSelect], hydro.θsMat[1:N_SoilSelect], hydro.ΨkgMac[1:N_SoilSelect], hydro.σMac[1:N_SoilSelect], θr_Psd[1:N_SoilSelect])

		else
			error("  $(option.psd.OptimizePsd) not found ")
		end # option.psd.OptimizePsd

		# STATISTICS NSE
		println("\n PSD MODEL:")
		Nse_Psd, Nse_Psd_Mean = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd[1:N_SoilSelect], Ψ_Rpart[1:N_SoilSelect,:], θ_Rpart[1:N_SoilSelect,:], hydro.θsMac[1:N_SoilSelect], hydro.θr[1:N_SoilSelect], hydro.ΨkgMat[1:N_SoilSelect], hydro.σMat[1:N_SoilSelect], hydro.θsMat[1:N_SoilSelect], hydro.ΨkgMac[1:N_SoilSelect], hydro.σMac[1:N_SoilSelect])
		# Nse_Psd_Mean = Statistics.mean(Nse_Psd[1:N_SoilSelect])


		# =================================================================================================================
		#       OPTIMIZATION HYDRAULIC PARAMETERS DERIVED FROM psd
		# =================================================================================================================
		if option.psd.HydroParam
		println("START Optimizing the hydraulic parameters derived from PSD")
		println("...")

		θsMat_Psd, θr_Psd_Kg, σMat_Psd, ΨkgMat_Psd, ~, θsMac_Psd, σMac_Psd, ΨkgMac_Psd, ~, Nse_θh_Uni_Psd, Nse_θh_Bim_Psd = optimiseCharacUnsat.OPTIMISE_CHARAC_UNSAT(Ψ_Rpart[1:N_SoilSelect,:], θ_Rpart[1:N_SoilSelect,:], N_Psd[1:N_SoilSelect], N_SoilSelect, hydro.θsMac[1:N_SoilSelect], θr_Psd=θr_Psd[1:N_SoilSelect]; Option_Data_Kθ = false )

		Nse_Psd_σ = 1. - stats.NASH_SUTCLIFE_MINIMIZE(hydro.σMat, σMat_Psd; Power=2.)
		Nse_Psd_Ψkg = 1. - stats.NASH_SUTCLIFE_MINIMIZE(hydro.ΨkgMat, ΨkgMat_Psd; Power=2.)
		ΔθsMacMat = hydro.θsMac[1:N_SoilSelect] .- hydro.θsMat[1:N_SoilSelect] 
		ΔθsMacMat_Psd = hydro.θsMac[1:N_SoilSelect] .- θsMat_Psd[1:N_SoilSelect] 
		Nse_Psd_θs = 1. - stats.NASH_SUTCLIFE_MINIMIZE(ΔθsMacMat, ΔθsMacMat_Psd; Power=2.)

		println("Nse_Psd_σ = $Nse_Psd_σ ; Nse_Psd_Ψkg = $Nse_Psd_Ψkg ; Nse_Psd_θs =$Nse_Psd_θs")

		table.SIMULATION(Nse_Psd_σ, Nse_Psd_Ψkg, Nse_Psd_θs, Nse_Psd_Mean, Nse_θΨ_Bim_Mean) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!mover a TABLE

		if  option.psd.Plot_σ_Ψkg # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!mover a PLOT
			plot.PLOTTING_σ_Ψkg(hydro.σMat[1:N_SoilSelect], σMat_Psd[1:N_SoilSelect], hydro.ΨkgMat[1:N_SoilSelect], ΨkgMat_Psd[1:N_SoilSelect], ΔθsMacMat[1:N_SoilSelect], ΔθsMacMat_Psd[1:N_SoilSelect])
		end # option.psd.Plot_σ_Ψkg

		println("END Optimizing the hydraulic parameters derived from PSD \n") 
		end # option.Psd.HydroParam =======================================================================================
		

		# =================================================================================================================
		#			TABLES
		# =================================================================================================================
		println("START WRITTING TABLE, \n")
		println(" ... ")

		if option.psd.OptimizePsd == "Single" 
			# table.SINGLEOPT_ξ1_ξ2(hydro.θsMac[1:N_SoilSelect], hydro.θr[1:N_SoilSelect], θr_Psd[1:N_SoilSelect], hydro.σMat[1:N_SoilSelect], hydro.ΨkgMat[1:N_SoilSelect], hydro.θsMat[1:N_SoilSelect], hydro.σMac[1:N_SoilSelect], hydro.ΨkgMac[1:N_SoilSelect], ξ1[1:N_SoilSelect], ξ2[1:N_SoilSelect], Nse_Psd[1:N_SoilSelect], Subclay[1:N_SoilSelect])
			table.SINGLEOPT_ξ1_ξ2(θr_Psd[1:N_SoilSelect], ξ1[1:N_SoilSelect], ξ2[1:N_SoilSelect], Nse_Psd[1:N_SoilSelect], Subclay[1:N_SoilSelect], hydro[1:N_SoilSelect]) #!!!!!!!!!!!!
		end

		if option.psd.HydroParam
			table.HYDRAULICPARAM_Psd(θsMat_Psd[1:N_SoilSelect], θr_Psd_Kg[1:N_SoilSelect], σMat_Psd[1:N_SoilSelect], ΨkgMat_Psd[1:N_SoilSelect], θsMac_Psd[1:N_SoilSelect], σMac_Psd[1:N_SoilSelect], ΨkgMac_Psd[1:N_SoilSelect], Nse_θh_Uni_Psd[1:N_SoilSelect],Nse_θh_Bim_Psd[1:N_SoilSelect], ∑Psd[1:N_SoilSelect,1:param.psd.N_Psd])
		end

		table.INTERGRANULARMIXING(hydro.θsMac[1:N_SoilSelect], hydro.θr[1:N_SoilSelect], θr_Psd[1:N_SoilSelect], hydro.σMat[1:N_SoilSelect], hydro.ΨkgMat[1:N_SoilSelect], hydro.θsMat[1:N_SoilSelect], hydro.σMac[1:N_SoilSelect], hydro.ΨkgMac[1:N_SoilSelect], ξ1[1:N_SoilSelect], ξ2[1:N_SoilSelect], Nse_Psd[1:N_SoilSelect], Subclay[1:N_SoilSelect], ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, 0.0, 0.0, 0.0)

		println("END WRITTING TABLE, \n")












				# ........................................................................


		# =================================================================================================================
		#       OPTIMIZATION SINGLE SOIL
		# =================================================================================================================
		function OPTIMIZATION_SINGLE_SOIL(iSoil, Psd, ∑Psd, Rpart, N_Psd, hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac, θr_Psd)
			if option.psd.SubclayOpt  && !(option.psd.∑Psd_2_ξ2) && !(option.psd.∑Psd_2_ξ1) # <><><><><><><>	

				Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(P[1], P[2], P[3], Rpart[1:N_Psd], N_Psd, Psd[1:N_Psd], ∑Psd[1:N_Psd], hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac, θr_Psd) ; SearchRange =[(param.psd.ξ1_Min, param.psd.ξ1_Max), (param.psd.ξ2_Min, param.psd.ξ2_Max), (param.psd.Wsubclay_Min, param.psd.Wsubclay_Max)], NumDimensions=3,  TraceMode=:silent)

				# Optimal values
				ξ1 = BlackBoxOptim.best_candidate(Optimization)[1]
				ξ2 = BlackBoxOptim.best_candidate(Optimization)[2]
				Subclay = BlackBoxOptim.best_candidate(Optimization)[3]
				Of_Psd = BlackBoxOptim.best_fitness(Optimization)

			elseif !(option.psd.SubclayOpt) && !(option.psd.∑Psd_2_ξ2) && !(option.psd.∑Psd_2_ξ1)# <><><><><><><># <><><><><><><>
				Subclay = param.psd.Subclay
		
				Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(P[1], P[2], Subclay, Rpart[1:N_Psd], N_Psd, Psd[1:N_Psd], ∑Psd[1:N_Psd], hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac, θr_Psd) ; SearchRange =[(param.psd.ξ1_Min, param.psd.ξ1_Max), (param.psd.ξ2_Min, param.psd.ξ2_Max)], NumDimensions=2,  TraceMode=:silent)

				# Optimal values
				ξ1 = BlackBoxOptim.best_candidate(Optimization)[1]
				ξ2 = BlackBoxOptim.best_candidate(Optimization)[2]
				Of_Psd = BlackBoxOptim.best_fitness(Optimization)

			elseif !(option.psd.SubclayOpt) && !(option.psd.∑Psd_2_ξ2)  && (option.psd.∑Psd_2_ξ1)
				Subclay = param.psd.Subclay
				ξ1 = param.psd.P_ξ1
		
				Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(ξ1 , P[1], Subclay, Rpart[1:N_Psd], N_Psd, Psd[1:N_Psd], ∑Psd[1:N_Psd], hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac, θr_Psd) ; SearchRange =[(param.psd.ξ2_Min, param.psd.ξ2_Max)], NumDimensions=1,  TraceMode=:silent)
				# Optimal values
				ξ2 = BlackBoxOptim.best_candidate(Optimization)[1]
				Of_Psd = BlackBoxOptim.best_fitness(Optimization)


			elseif option.psd.∑Psd_2_ξ2  && !(option.psd.∑Psd_2_ξ1)
				Subclay = param.psd.Subclay
				ξ2 = psdFunc.∑PSD_2_ξ2(∑Psd[param.psd.∑Psd_2_ξ2_Size])

				Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(P[1], ξ2 , Subclay, Rpart[1:N_Psd], N_Psd, Psd[1:N_Psd], ∑Psd[1:N_Psd], hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac, θr_Psd) ; SearchRange =[(param.psd.ξ1_Min, param.psd.ξ1_Max)], NumDimensions=1,  TraceMode=:silent)

				# Optimal values
				ξ1 = BlackBoxOptim.best_candidate(Optimization)[1]
				Of_Psd = BlackBoxOptim.best_fitness(Optimization)


			elseif option.psd.∑Psd_2_ξ2  && option.psd.∑Psd_2_ξ1
				Subclay = param.psd.Subclay # Average value

				ξ2 =  psdFunc.∑PSD_2_ξ2(∑Psd[param.psd.∑Psd_2_ξ2_Size])
				ξ1 = param.psd.P_ξ1
				Of_Psd = OF_SINGLE_SOIL(ξ1, ξ2, param.psd.Subclay, Rpart[1:N_Psd], N_Psd, Psd[1:N_Psd], ∑Psd[1:N_Psd], hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac, θr_Psd)
			else
				error(" NO SUITABLE OPTIONS FOUND ")
			end # if option.psd.

			# Recording the optimal parameters
			θ_Rpart, Ψ_Rpart = psdFunc.PSD_MODEL(Rpart, Psd, ∑Psd, N_Psd, hydro.θsMac, θr_Psd, Subclay, ξ1, ξ2)

			# Correction of Subclay  and deriving the corrected Psd and  ∑Psd
			Psd, ∑Psd = psdFunc.SUBCLAY_CORRECTION(∑Psd, Subclay, N_Psd)

			return ξ1, ξ2, Ψ_Rpart, θ_Rpart, Of_Psd, Psd, ∑Psd, Subclay
		end # OPTIMIZATION_SINGLE_SOIL=====================================================================================


		function OF_SINGLE_SOIL(ξ1, ξ2, Subclay, Rpart, N_Psd, Psd, ∑Psd, hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac, θr_Psd) 
			θ_Rpart = zeros(Float64, N_Psd)
			Ψ_Rpart = zeros(Float64, N_Psd)
			θΨ = zeros(Float64, N_Psd)

			θ_Rpart[1:N_Psd], Ψ_Rpart[1:N_Psd] = psdFunc.PSD_MODEL(Rpart, Psd, ∑Psd, N_Psd, hydro.θsMac, θr_Psd, Subclay, ξ1, ξ2)

			# For every class
			Of = 0.0
			for iRpart = 1:N_Psd
				# Observed data
				θΨ[iRpart] = wrc.kg.Ψ_2_θdual(Ψ_Rpart[iRpart], hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac)
				# Summing the error
				Of += (θΨ[iRpart] - θ_Rpart[iRpart])^2.0
			end
			Of = (Of / N_Psd) ^ 0.5
			return Of
		end # OF_SINGLE_SOIL ===============

		
		# =================================================================================================================
		#	Optimisation ALL SOILS
		# ================================================================================================================= 
		function OPTIMIZATION_ALL_SOIL(N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac, θr_Psd)

			function OF_ALL_SOIL(ξ1, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac, θr_Psd)
				Of_AllSoil = 0.
				@simd for iSoil = 1:N_SoilSelect
					ξ2 =  psdFunc.∑PSD_2_ξ2(∑Psd[iSoil, param.psd.∑Psd_2_ξ2_Size]; ∑Psd_2_ξ2_β1=∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=∑Psd_2_ξ2_β2) #####################for new table model 4######
					# ξ2 = ∑Psd_2_ξ2_β1   ############################for new table model 1######

					Of_AllSoil += OF_SINGLE_SOIL(ξ1, ξ2, Subclay, Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], hydro.θsMac[iSoil], hydro.θr[iSoil], hydro.ΨkgMat[iSoil], hydro.σMat[iSoil], hydro.θsMat[iSoil], hydro.ΨkgMac[iSoil], hydro.σMac[iSoil], θr_Psd[iSoil])
				end # for loop
				return Of_AllSoil
			end # function OF_ALL_SOILS


			if option.psd.SubclayOpt 
				Optimization = BlackBoxOptim.bboptimize(P -> OF_ALL_SOIL(P[1], P[2], P[3], P[4], N_SoilSelect, Psd[1:N_SoilSelect,:], ∑Psd[1:N_SoilSelect,:], Rpart[1:N_SoilSelect,:], N_Psd[1:N_SoilSelect], hydro.θsMac[1:N_SoilSelect], hydro.θr[1:N_SoilSelect], hydro.ΨkgMat[1:N_SoilSelect], hydro.σMat[1:N_SoilSelect], hydro.θsMat[1:N_SoilSelect], hydro.ΨkgMac[1:N_SoilSelect], hydro.σMac[1:N_SoilSelect], θr_Psd[1:N_SoilSelect]);  SearchRange =[(param.psd.ξ1_Min, param.psd.ξ1_Max), (param.psd.∑Psd_2_ξ2_β1_Min,param.psd.∑Psd_2_ξ2_β1_Max), (param.psd.∑Psd_2_ξ2_β2_Min,param.psd.∑Psd_2_ξ2_β2_Max), (param.psd.Wsubclay_Min, param.psd.Wsubclay_Max)], NumDimensions=4, TraceMode=:silent )

				ξ1_All =  BlackBoxOptim.best_candidate(Optimization)[1]
				∑Psd_2_ξ2_β1 = BlackBoxOptim.best_candidate(Optimization)[2]
				∑Psd_2_ξ2_β2 =  BlackBoxOptim.best_candidate(Optimization)[3]
				Wsubclay_All =  BlackBoxOptim.best_candidate(Optimization)[4]

			elseif !(option.psd.SubclayOpt)
				Wsubclay_All = param.psd.Subclay

				Optimization = BlackBoxOptim.bboptimize(P -> OF_ALL_SOIL(P[1], P[2], P[3], Wsubclay_All, N_SoilSelect, Psd[1:N_SoilSelect,:], ∑Psd[1:N_SoilSelect,:], Rpart[1:N_SoilSelect,:], N_Psd[1:N_SoilSelect], hydro.θsMac[1:N_SoilSelect], hydro.θr[1:N_SoilSelect], hydro.ΨkgMat[1:N_SoilSelect], hydro.σMat[1:N_SoilSelect], hydro.θsMat[1:N_SoilSelect], hydro.ΨkgMac[1:N_SoilSelect], hydro.σMac[1:N_SoilSelect], θr_Psd[1:N_SoilSelect]);  SearchRange =[(param.psd.ξ1_Min, param.psd.ξ1_Max), (param.psd.∑Psd_2_ξ2_β1_Min,param.psd.∑Psd_2_ξ2_β1_Max), (param.psd.∑Psd_2_ξ2_β2_Min,param.psd.∑Psd_2_ξ2_β2_Max)], NumDimensions=3, TraceMode=:silent )

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

				θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] = psdFunc.PSD_MODEL(Rpart[iSoil,1:N_Psd[iSoil]], Psd[iSoil,:], ∑Psd[iSoil,:], N_Psd[iSoil], hydro.θsMac[iSoil], θr_Psd[iSoil], Subclay[iSoil], ξ1[iSoil], ξ2[iSoil])

				Of_Psd[iSoil] = OF_SINGLE_SOIL(ξ1[iSoil], ξ2[iSoil], Subclay[iSoil], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], hydro.θsMac[iSoil], hydro.θr[iSoil], hydro.ΨkgMat[iSoil], hydro.σMat[iSoil], hydro.θsMat[iSoil], hydro.ΨkgMac[iSoil], hydro.σMac[iSoil], θr_Psd[iSoil])
			end # For

			return ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd
		end # function OPTIMIZATION_ALL_SOIL ==============================================================================

		
		# ==================================================================================================================
		#	Parameters ALL SOILS
		# ==================================================================================================================
		function PARAMETERS_ALL_SOIL(N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac, θr_Psd)

			function OF_ALL_SOIL(ξ1, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, hydro.θsMac, hydro.θr, hydro.ΨkgMat, hydro.σMat, hydro.θsMat, hydro.ΨkgMac, hydro.σMac, θr_Psd)
				Of_AllSoil = 0.0
				@simd for iSoil = 1:N_SoilSelect
					ξ2 =  psdFunc.∑PSD_2_ξ2(∑Psd[iSoil, param.psd.∑Psd_2_ξ2_Size]; ∑Psd_2_ξ2_β1=param.psd.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=param.psd.∑Psd_2_ξ2_β2) 				

					Of_AllSoil += OF_SINGLE_SOIL(ξ1, ξ2, Subclay, Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], hydro.θsMac[iSoil], hydro.θr[iSoil], hydro.ΨkgMat[iSoil], hydro.σMat[iSoil], hydro.θsMat[iSoil], hydro.ΨkgMac[iSoil], hydro.σMac[iSoil], θr_Psd[iSoil])
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

				θ_Rpart[iSoil,1:N_Psd[iSoil]], Ψ_Rpart[iSoil,1:N_Psd[iSoil]] = psdFunc.PSD_MODEL(Rpart[iSoil,1:N_Psd[iSoil]], Psd[iSoil,:], ∑Psd[iSoil,:], N_Psd[iSoil], hydro.θsMac[iSoil], θr_Psd[iSoil], Subclay[iSoil], ξ1[iSoil], ξ2[iSoil])

				Of_Psd[iSoil] = OF_SINGLE_SOIL(ξ1[iSoil], ξ2[iSoil], Subclay[iSoil], Rpart[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil], Psd[iSoil,1:N_Psd[iSoil]], ∑Psd[iSoil,1:N_Psd[iSoil]], hydro.θsMac[iSoil], hydro.θr[iSoil], hydro.ΨkgMat[iSoil], hydro.σMat[iSoil], hydro.θsMat[iSoil], hydro.ΨkgMac[iSoil], hydro.σMac[iSoil], θr_Psd[iSoil])
			end # For

			return ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd
		end # function PARAMETERS_ALL_SOIL =================================================================================

		

	end # function PSD_MAIN


end # module PSD