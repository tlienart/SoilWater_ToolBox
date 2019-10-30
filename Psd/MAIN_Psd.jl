module psd
import ..option, ..param, ..wrc, ..kunsat, ..cst, ..path, ..stats, ..psdFunc



	using Statistics, BlackBoxOptim, JuliaDB

	# ======================================================================================
	#          PSD_START
	# ======================================================================================
	function PSD_START(Nsample, θsMat, θr, σMat, ΨkgMat, KsMat, θsMac, σMac, ΨkgMac, KsMac, Ψ_θΨ, θ_θΨ, N_θΨ, Nse_θΨ_Bim_Mean, Flag_Good ; N_Kθ=ones(Int8, Nsample), K_Kθ=zeros(Float64,Nsample,1), Ψ_Kθ=zeros(Float64,Nsample,1))

		println("option.SubclayOpt = $(option.SubclayOpt), option.∑Psd_2_ξ2 = $(option.∑Psd_2_ξ2), option.∑Psd_2_ξ1 = $(option.∑Psd_2_ξ1), \n")

		plot.PLOTTING_MIXING_PARAM()
		plot.PLOTTING_THETA_R_MODEL()
		plot.PLOTTING_R_2_PSI_MODEL()
		# plot.PLOTTING_ξ2_MODEL()


		# # Number of soils wanted to run
		# 	Nsample = min(Nsample, param.N_Soil_Select)

		# RESERVING MEMORY
			local Npsd = 100 # Just a guess
			Ψ_Rpart = zeros(Float64, (Nsample,Npsd))
			ξ2 = zeros(Float64, Nsample)
			ξ1 = zeros(Float64, Nsample)
			ξ = zeros(Float64, Npsd) # Problem need to solve
			θr_Psd =  zeros(Float64,Nsample)
			θ_Rpart = zeros(Float64, (Nsample,Npsd))
			Subclay = zeros(Float64, Nsample)
			Rpart_Whole = zeros(Float64, Npsd)
			Rpart = zeros(Float64, (Nsample,Npsd))
			Psd = zeros(Float64, (Nsample,Npsd))
			Of_Psd = zeros(Float64, Nsample)
			Nrpart = zeros(Int32, Nsample)
			∑Psd = zeros(Float64, (Nsample,Npsd))
			Nse_Psd = zeros(Float64, Nsample)
			Nse_Psd_θh = zeros(Float64, Nsample)
			θr_Psd_Kg = zeros(Float64, Nsample)

		# INITIAL VALUES WHICH WILL CHANGE
			∑Psd_2_ξ2_β1 = param.∑Psd_2_ξ2_β1
			∑Psd_2_ξ2_β2 = param.∑Psd_2_ξ2_β2

		# Computing the radius of the particles from the diameters
			Rpart_Whole = param.Dpart ./ 2.0

		# Reading PSD data
			∑Psd = reading.PSD(path.Psd, Flag_Good)

		# =================== FOR EVERY SOIL SAMPLES ===================
		for iSoil=1:Nsample
			# Compute maximum N particle size for each iSoil
			# Write Rpart such that it varies for every iSoil
			for iPsd = 1:length(Rpart)
					Rpart[iSoil,iPsd] = Rpart_Whole[iPsd]
					if ∑Psd[iSoil, iPsd] >= 0.99999
						Nrpart_Count = iPsd
						Nrpart[iSoil] = Nrpart_Count
						break
					end
				end
			
				
			# Compute PSD from ∑PSD
			Psd[iSoil,1:Nrpart[iSoil]] = psdFunc.∑PSD_2_PSD(∑Psd[iSoil,1:Nrpart[iSoil]], Nrpart[iSoil])

		end # looping over soils


		# =================== COMPUTTING  θr model from Psd ===================
			if option.Psd_2_θr == "Opt"
				θr_Psd = psdFunc.OPTIMIZE_PSD_2_θr(Nsample, θr[1:Nsample], ∑Psd[1:Nsample,:])

			elseif  option.Psd_2_θr == "Cst"
				@simd for iSoil=1:Nsample
					θr_Psd[iSoil] = param.θr_Cst
				end
				
			elseif  option.Psd_2_θr == "Param"
				println("Optimize θr = Psd_2_θr_α1 = $(param.Psd_2_θr_α1) ; param.Psd_2_θr_α2 = $(param.Psd_2_θr_α2)")
				@simd for iSoil=1:Nsample
					θr_Psd[iSoil] = psdFunc.PSD_2_θr(∑Psd[iSoil,param.Psd_2_θr_Size])
				end
			else
				error("option.Psd_2_θr = $option.Psd_2_θr  not allowed option.Psd_2_θr must be either (1)'Opt' or (2) 'Cst' or (3) 'Param' ")
			end # if option.Psd_2_θr

			Nse_θr = 1.0 - stats.NASH_SUTCLIFE_MINIMIZE(θr[1:Nsample], θr_Psd[1:Nsample])

			println("NASH_SUTCLIFFE θr_Psd = $Nse_θr \n")

			if option.Plot_θr
				plot.PLOTTING_θr(θr[1:Nsample], θr_Psd[1:Nsample], ∑Psd[1:Nsample,param.Psd_2_θr_Size])
			end # option.Plot_θr

		# ........................................................................


		if option.OptimizePsd == "Single" # <> <> <> <> <> <>
			∑Of_Psd = 0.0

			@simd for iSoil=1:Nsample
				ξ1[iSoil], ξ2[iSoil], Ψ_Rpart[iSoil,1:Nrpart[iSoil]], θ_Rpart[iSoil,1:Nrpart[iSoil]], Of_Psd[iSoil], Psd[iSoil,1:Nrpart[iSoil]], ∑Psd[iSoil,1:Nrpart[iSoil]], Subclay[iSoil] = OPTIMIZATION_SINGLE_SOIL(iSoil, Psd[iSoil,1:Nrpart[iSoil]], ∑Psd[iSoil,1:Nrpart[iSoil]], Rpart[iSoil,1:Nrpart[iSoil]], Nrpart[iSoil], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil], θr_Psd[iSoil])
			end # Loop single

		elseif option.OptimizePsd == "All" # <> <> <> <> <> <>
			∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd = OPTIMIZATION_ALL_SOIL(Nsample, Psd[1:Nsample,:], ∑Psd[1:Nsample,:], Rpart[1:Nsample,:], Nrpart[1:Nsample], θsMac[1:Nsample], θr[1:Nsample], ΨkgMat[1:Nsample], σMat[1:Nsample], θsMat[1:Nsample], ΨkgMac[1:Nsample], σMac[1:Nsample], θr_Psd[1:Nsample])
		
		elseif option.OptimizePsd == "Run" # <> <> <> <> <> <>
			∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd = PARAMETERS_ALL_SOIL(Nsample, Psd[1:Nsample,:], ∑Psd[1:Nsample,:], Rpart[1:Nsample,:], Nrpart[1:Nsample], θsMac[1:Nsample], θr[1:Nsample], ΨkgMat[1:Nsample], σMat[1:Nsample], θsMat[1:Nsample], ΨkgMac[1:Nsample], σMac[1:Nsample], θr_Psd[1:Nsample])

		else
			error("  $(option.OptimizePsd) not found ")
		end # option.OptimizePsd

		# STATISTICS NSE
			println("\n PSD MODEL:")
			Nse_Psd, Nse_Psd_Mean = stats.NASH_SUTCLIFFE_θΨ(Nsample, Nrpart[1:Nsample], Ψ_Rpart[1:Nsample,:], θ_Rpart[1:Nsample,:], θsMac[1:Nsample], θr[1:Nsample], ΨkgMat[1:Nsample], σMat[1:Nsample], θsMat[1:Nsample], ΨkgMac[1:Nsample], σMac[1:Nsample])

			# Nse_Psd_Mean = Statistics.mean(Nse_Psd[1:Nsample])


		# =================================================================================================================
		#       OPTIMIZATION HYDRAULIC PARAMETERS DERIVED FROM psd
		# =================================================================================================================
		if option.Psd.HydrauParam
			println("START Optimizing the hydraulic parameters derived from PSD")
			println("...")

			θsMat_Psd, θr_Psd_Kg, σMat_Psd, ΨkgMat_Psd, ~, θsMac_Psd, σMac_Psd, ΨkgMac_Psd, ~, Nse_θh_Uni_Psd, Nse_θh_Bim_Psd = optimiseCharacUnsat.OPTIMISE_CHARAC_UNSAT(Ψ_Rpart[1:Nsample,:], θ_Rpart[1:Nsample,:], Nrpart[1:Nsample], Nsample, θsMac[1:Nsample], θr_Psd=θr_Psd[1:Nsample]; Option_Data_Kθ = false )

			Nse_Psd_σ = 1. - stats.NASH_SUTCLIFE_MINIMIZE(σMat, σMat_Psd; Power=2.)
			Nse_Psd_Ψkg = 1. - stats.NASH_SUTCLIFE_MINIMIZE(ΨkgMat, ΨkgMat_Psd; Power=2.)
			ΔθsMacMat = θsMac[1:Nsample] .- θsMat[1:Nsample] 
			ΔθsMacMat_Psd = θsMac[1:Nsample] .- θsMat_Psd[1:Nsample] 
			Nse_Psd_θs = 1. - stats.NASH_SUTCLIFE_MINIMIZE(ΔθsMacMat, ΔθsMacMat_Psd ; Power=2.)

			println("Nse_Psd_σ = $Nse_Psd_σ ; Nse_Psd_Ψkg = $Nse_Psd_Ψkg ; Nse_Psd_θs =$Nse_Psd_θs")

			table.SIMULATION(Nse_Psd_σ, Nse_Psd_Ψkg, Nse_Psd_θs, Nse_Psd_Mean, Nse_θΨ_Bim_Mean)

			if  option.Plot_σ_Ψkg
				plot.PLOTTING_σ_Ψkg(σMat[1:Nsample], σMat_Psd[1:Nsample], ΨkgMat[1:Nsample], ΨkgMat_Psd[1:Nsample], ΔθsMacMat[1:Nsample], ΔθsMacMat_Psd[1:Nsample])
			end # option.Plot_σ_Ψkg

			println("END Optimizing the hydraulic parameters derived from PSD \n") 
		end # option.Psd.HydrauParam
	

		# =================================================================================================================
		#			TABLES
		# =================================================================================================================
			println("START WRITTING TABLE, \n")
			println(" ... ")

			if option.OptimizePsd == "Single" 
				table.SINGLEOPT_T1_T2(θsMac[1:Nsample], θr[1:Nsample], θr_Psd[1:Nsample], σMat[1:Nsample], ΨkgMat[1:Nsample], θsMat[1:Nsample], σMac[1:Nsample], ΨkgMac[1:Nsample], ξ1[1:Nsample], ξ2[1:Nsample], Nse_Psd[1:Nsample], Subclay[1:Nsample])
			end

			if option.Psd_2_HydrauParam
				table.HYDRAULICPARAM_Psd(θsMat_Psd[1:Nsample], θr_Psd_Kg[1:Nsample], σMat_Psd[1:Nsample], ΨkgMat_Psd[1:Nsample], θsMac_Psd[1:Nsample], σMac_Psd[1:Nsample], ΨkgMac_Psd[1:Nsample], Nse_θh_Uni_Psd[1:Nsample],Nse_θh_Bim_Psd[1:Nsample], ∑Psd[1:Nsample,1:param.N_Psd])
			end

			table.INTERGRANULARMIXING(θsMac[1:Nsample], θr[1:Nsample], θr_Psd[1:Nsample], σMat[1:Nsample], ΨkgMat[1:Nsample], θsMat[1:Nsample], σMac[1:Nsample], ΨkgMac[1:Nsample], ξ1[1:Nsample], ξ2[1:Nsample], Nse_Psd[1:Nsample], Subclay[1:Nsample], ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, 0.0, 0.0, 0.0)

			println("END WRITTING TABLE, \n")


		# =================================================================================================================
		#			PLOTTING
		# =================================================================================================================
			
		plot.PLOTTING_MIXING_PARAM(;ξ1=6.56, ξ2=0.18)
		
			if option.Plot_∑Psd_2_ξ2
				Data = JuliaDB.loadtable(path.Table * "\\IndividualSamples_T1_T2.csv")
				ξ2_Obs = JuliaDB.select(Data,:Tau2)
				plot.PLOTTING_∑Psd_2_ξ2(ξ2_Obs[1:Nsample], ξ2[1:Nsample], ∑Psd[1:Nsample, param.∑Psd_2_ξ2_Size])

				#  Paper graph
				plot.PLOTTING_ξ2_MODEL(ξ2_Obs[1:Nsample], ξ2[1:Nsample], ∑Psd[1:Nsample, param.∑Psd_2_ξ2_Size])

				# Statistics
				Nse_ξ2 = 1.0 - stats.NASH_SUTCLIFE_MINIMIZE(ξ2_Obs[1:Nsample], ξ2[1:Nsample])
				println("NSE_ξ2 = $Nse_ξ2")

			end # option.Plot_∑Psd_2_ξ2

			if option.Plot_OFsingle_OFmeas
				Data = JuliaDB.loadtable(path.Table * "\\IndividualSamples_T1_T2.csv")
				Nse_T1T2 = JuliaDB.select(Data,:NSE_Psd)

				Data = JuliaDB.loadtable(path.OptimizeKΨ)
				Nse_Obs = JuliaDB.select(Data,:Nse_THETAh_Bim)

				plot.PLOTTING_OF(Nse_T1T2[1:Nsample], Nse_Obs[1:Nsample], Nse_Psd[1:Nsample])
			end # option.Plot_OFsingle_OFmeas

			if option.Plotting
				@simd for iSoil = 1:Nsample
					println("Plotting soil $iSoil")

					#   INTERGRANULARMIXING
					@simd for iRpart = 1: Nrpart[iSoil]
						ξ[iRpart] = psdFunc.INTERGRANULARMIXING(Rpart[iSoil,iRpart], ξ1[iSoil], ξ2[iSoil])
					end

					if option.OptimizeKΨ && option.Psd_2_HydrauParam
						plot.PLOTTING(iSoil, Nrpart[iSoil], Rpart[iSoil,1:Nrpart[iSoil]], Psd[iSoil,1:Nrpart[iSoil]], ∑Psd[iSoil,1:Nrpart[iSoil]], θsMac[iSoil], θr[iSoil], σMat[iSoil], ΨkgMat[iSoil], KsMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMac[iSoil], ΨkgMac[iSoil], Ψ_θΨ[iSoil,1:N_θΨ[iSoil]], θ_θΨ[iSoil,1:N_θΨ[iSoil]], N_θΨ[iSoil], Ψ_Rpart[iSoil,1:Nrpart[iSoil]], θ_Rpart[iSoil,1:Nrpart[iSoil]], ξ[1:Nrpart[iSoil]]; N_Kθ=N_Kθ[iSoil], K_Kθ=K_Kθ[iSoil,1:N_Kθ[iSoil]], Ψ_Kθ=Ψ_Kθ[iSoil,1:N_Kθ[iSoil]],θsMat_Psd=θsMat_Psd[iSoil], θr_Psd=θr_Psd[iSoil], σMat_Psd=σMat_Psd[iSoil], ΨkgMat_Psd=ΨkgMat_Psd[iSoil], θsMac_Psd=θsMac_Psd[iSoil], σMac_Psd=σMac_Psd[iSoil], ΨkgMac_Psd=ΨkgMac_Psd[iSoil])

					elseif option.OptimizeKΨ && !option.Psd_2_HydrauParam
						plot.PLOTTING(iSoil, Nrpart[iSoil], Rpart[iSoil,1:Nrpart[iSoil]], Psd[iSoil,1:Nrpart[iSoil]], ∑Psd[iSoil,1:Nrpart[iSoil]], θsMac[iSoil], θr[iSoil], σMat[iSoil], ΨkgMat[iSoil], KsMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMac[iSoil], ΨkgMac[iSoil], Ψ_θΨ[iSoil,1:N_θΨ[iSoil]], θ_θΨ[iSoil,1:N_θΨ[iSoil]], N_θΨ[iSoil], Ψ_Rpart[iSoil,1:Nrpart[iSoil]], θ_Rpart[iSoil,1:Nrpart[iSoil]], ξ[1:Nrpart[iSoil]]; N_Kθ=N_Kθ[iSoil], K_Kθ=K_Kθ[iSoil,1:N_Kθ[iSoil]], Ψ_Kθ=Ψ_Kθ[iSoil,1:N_Kθ[iSoil]])

					elseif !option.OptimizeKΨ && option.Psd_2_HydrauParam
						plot.PLOTTING(iSoil, Nrpart[iSoil], Rpart[iSoil,1:Nrpart[iSoil]],Psd[iSoil,1:Nrpart[iSoil]], ∑Psd[iSoil,1:Nrpart[iSoil]], θsMac[iSoil], θr[iSoil], σMat[iSoil], ΨkgMat[iSoil], KsMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMac[iSoil], ΨkgMac[iSoil], Ψ_θΨ[iSoil,1:N_θΨ[iSoil]], θ_θΨ[iSoil,1:N_θΨ[iSoil]], N_θΨ[iSoil], Ψ_Rpart[iSoil,1:Nrpart[iSoil]], θ_Rpart[iSoil,1:Nrpart[iSoil]], ξ[1:Nrpart[iSoil]];θsMat_Psd=θsMat_Psd[iSoil], θr_Psd=θr_Psd[iSoil], σMat_Psd=σMat_Psd[iSoil], ΨkgMat_Psd=ΨkgMat_Psd[iSoil], θsMac_Psd=θsMac_Psd[iSoil], σMac_Psd=σMac_Psd[iSoil], ΨkgMac_Psd=ΨkgMac_Psd[iSoil])
					
					elseif !option.OptimizeKΨ && !option.Psd_2_HydrauParam
						plot.PLOTTING(iSoil, Nrpart[iSoil], Rpart[iSoil,1:Nrpart[iSoil]], Psd[iSoil,1:Nrpart[iSoil]], ∑Psd[iSoil,1:Nrpart[iSoil]], θsMac[iSoil], θr[iSoil], σMat[iSoil], ΨkgMat[iSoil], KsMat[iSoil], KsMac[iSoil], θsMat[iSoil], σMac[iSoil], ΨkgMac[iSoil], Ψ_θΨ[iSoil,1:N_θΨ[iSoil]], θ_θΨ[iSoil,1:N_θΨ[iSoil]], N_θΨ[iSoil], Ψ_Rpart[iSoil,1:Nrpart[iSoil]], θ_Rpart[iSoil,1:Nrpart[iSoil]], ξ[1:Nrpart[iSoil]])
					end # elseif

				end # looping over soils
			end #  if option.Plotting
		return
	end # function PSD_START



	# =================================================================================================================
	#       OPTIMIZATION SINGLE SOIL
	# =================================================================================================================
	function OPTIMIZATION_SINGLE_SOIL(iSoil, Psd, ∑Psd, Rpart, Nrpart, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)
		if option.SubclayOpt  && !(option.∑Psd_2_ξ2) && !(option.∑Psd_2_ξ1) # <><><><><><><>	

			Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(P[1], P[2], P[3], Rpart[1:Nrpart], Nrpart, Psd[1:Nrpart], ∑Psd[1:Nrpart], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd) ; SearchRange =[(param.ξ1_Min, param.ξ1_Max), (param.ξ2_Min, param.ξ2_Max), (param.Wsubclay_Min, param.Wsubclay_Max)], NumDimensions=3,  TraceMode=:silent)

			 # Optimal values
			ξ1 = BlackBoxOptim.best_candidate(Optimization)[1]
			ξ2 = BlackBoxOptim.best_candidate(Optimization)[2]
			Subclay = BlackBoxOptim.best_candidate(Optimization)[3]
			Of_Psd = BlackBoxOptim.best_fitness(Optimization)

		elseif !(option.SubclayOpt) && !(option.∑Psd_2_ξ2) && !(option.∑Psd_2_ξ1)# <><><><><><><># <><><><><><><>
			Subclay = param.Subclay
	
			Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(P[1], P[2], Subclay , Rpart[1:Nrpart], Nrpart, Psd[1:Nrpart], ∑Psd[1:Nrpart], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd) ; SearchRange =[(param.ξ1_Min, param.ξ1_Max), (param.ξ2_Min, param.ξ2_Max)], NumDimensions=2,  TraceMode=:silent)

			 # Optimal values
			ξ1 = BlackBoxOptim.best_candidate(Optimization)[1]
			ξ2 = BlackBoxOptim.best_candidate(Optimization)[2]
			Of_Psd = BlackBoxOptim.best_fitness(Optimization)

		elseif !(option.SubclayOpt) && !(option.∑Psd_2_ξ2)  && (option.∑Psd_2_ξ1)
			Subclay = param.Subclay

			ξ1 = param.P_ξ1
	
			Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(ξ1 , P[1], Subclay , Rpart[1:Nrpart], Nrpart, Psd[1:Nrpart], ∑Psd[1:Nrpart], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd) ; SearchRange =[(param.ξ2_Min, param.ξ2_Max)], NumDimensions=1,  TraceMode=:silent)
			 # Optimal values
			ξ2 = BlackBoxOptim.best_candidate(Optimization)[1]
			Of_Psd = BlackBoxOptim.best_fitness(Optimization)


		elseif option.∑Psd_2_ξ2  && !(option.∑Psd_2_ξ1)
			Subclay = param.Subclay
			ξ2 = psdFunc.∑PSD_2_ξ2(∑Psd[param.∑Psd_2_ξ2_Size])

			Optimization = BlackBoxOptim.bboptimize(P -> OF_SINGLE_SOIL(P[1], ξ2 , Subclay , Rpart[1:Nrpart], Nrpart, Psd[1:Nrpart], ∑Psd[1:Nrpart], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd) ; SearchRange =[(param.ξ1_Min, param.ξ1_Max)], NumDimensions=1,  TraceMode=:silent)

			 # Optimal values
			ξ1 = BlackBoxOptim.best_candidate(Optimization)[1]
			Of_Psd = BlackBoxOptim.best_fitness(Optimization)


		elseif option.∑Psd_2_ξ2  && option.∑Psd_2_ξ1
			Subclay = param.Subclay # Average value

			ξ2 =  psdFunc.∑PSD_2_ξ2(∑Psd[param.∑Psd_2_ξ2_Size])

			ξ1 = param.P_ξ1

			Of_Psd = OF_SINGLE_SOIL(ξ1, ξ2, param.Subclay, Rpart[1:Nrpart], Nrpart, Psd[1:Nrpart], ∑Psd[1:Nrpart], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)
		else
			error(" NO SUITABLE OPTIONS FOUND ")
		end # if option.

		# Recording the optimal parameters
		θ_Rpart, Ψ_Rpart = psdFunc.PSD_MODEL(Rpart, Psd, ∑Psd, Nrpart, θsMac, θr_Psd, Subclay, ξ1, ξ2)

		# Correction of Subclay  and deriving the corrected Psd and  ∑Psd
		Psd, ∑Psd = psdFunc.SUBCLAY_CORRECTION(∑Psd, Subclay, Nrpart)

		return ξ1, ξ2, Ψ_Rpart, θ_Rpart, Of_Psd, Psd, ∑Psd, Subclay
	end # OPTIMIZATION_SINGLE_SOIL



	function OF_SINGLE_SOIL(ξ1, ξ2, Subclay, Rpart, Nrpart, Psd, ∑Psd, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd) 
		θ_Rpart = zeros(Float64, Nrpart)
		Ψ_Rpart = zeros(Float64, Nrpart)
		θΨ = zeros(Float64, Nrpart)

		θ_Rpart[1:Nrpart], Ψ_Rpart[1:Nrpart] = psdFunc.PSD_MODEL(Rpart, Psd, ∑Psd, Nrpart, θsMac, θr_Psd, Subclay, ξ1, ξ2)

		# For every class
		Of = 0.0
		for iRpart = 1:Nrpart
			# Observed data
			θΨ[iRpart] = wrc.kg.Ψ_2_θdual(Ψ_Rpart[iRpart], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac)
			# Summing the error
			Of += (θΨ[iRpart] - θ_Rpart[iRpart])^2.0
		end
		Of = (Of / Nrpart) ^ 0.5
		return Of
	end # OF_SINGLE_SOIL ===============



	#= =================================================================================================================
	       Optimisation ALL SOILS
	================================================================================================================= =#
	 function OPTIMIZATION_ALL_SOIL(Nsample, Psd, ∑Psd, Rpart, Nrpart, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)

		function OF_ALL_SOIL(ξ1, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, Nsample, Psd, ∑Psd, Rpart, Nrpart, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)
			Of_AllSoil = 0.
			@simd for iSoil = 1:Nsample
				ξ2 =  psdFunc.∑PSD_2_ξ2(∑Psd[iSoil, param.∑Psd_2_ξ2_Size]; ∑Psd_2_ξ2_β1=∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=∑Psd_2_ξ2_β2) #####################for new table model 4######
				# ξ2 = ∑Psd_2_ξ2_β1   ############################for new table model 1######

				Of_AllSoil += OF_SINGLE_SOIL(ξ1, ξ2, Subclay, Rpart[iSoil,1:Nrpart[iSoil]], Nrpart[iSoil], Psd[iSoil,1:Nrpart[iSoil]], ∑Psd[iSoil,1:Nrpart[iSoil]], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil], θr_Psd[iSoil])
			end # for loop
			return Of_AllSoil
		end # function OF_ALL_SOILS


		if option.SubclayOpt 
			Optimization = BlackBoxOptim.bboptimize(P -> OF_ALL_SOIL(P[1], P[2], P[3], P[4], Nsample, Psd[1:Nsample,:], ∑Psd[1:Nsample,:], Rpart[1:Nsample,:], Nrpart[1:Nsample], θsMac[1:Nsample], θr[1:Nsample], ΨkgMat[1:Nsample], σMat[1:Nsample], θsMat[1:Nsample], ΨkgMac[1:Nsample], σMac[1:Nsample], θr_Psd[1:Nsample]);  SearchRange =[(param.ξ1_Min, param.ξ1_Max), (param.∑Psd_2_ξ2_β1_Min,param.∑Psd_2_ξ2_β1_Max), (param.∑Psd_2_ξ2_β2_Min,param.∑Psd_2_ξ2_β2_Max), (param.Wsubclay_Min, param.Wsubclay_Max)], NumDimensions=4, TraceMode=:silent )

			ξ1_All =  BlackBoxOptim.best_candidate(Optimization)[1]
			∑Psd_2_ξ2_β1 = BlackBoxOptim.best_candidate(Optimization)[2]
			∑Psd_2_ξ2_β2 =  BlackBoxOptim.best_candidate(Optimization)[3]
			Wsubclay_All =  BlackBoxOptim.best_candidate(Optimization)[4]

		elseif !(option.SubclayOpt)
			Wsubclay_All = param.Subclay

			Optimization = BlackBoxOptim.bboptimize(P -> OF_ALL_SOIL(P[1], P[2], P[3], Wsubclay_All, Nsample, Psd[1:Nsample,:], ∑Psd[1:Nsample,:], Rpart[1:Nsample,:], Nrpart[1:Nsample], θsMac[1:Nsample], θr[1:Nsample], ΨkgMat[1:Nsample], σMat[1:Nsample], θsMat[1:Nsample], ΨkgMac[1:Nsample], σMac[1:Nsample], θr_Psd[1:Nsample]);  SearchRange =[(param.ξ1_Min, param.ξ1_Max), (param.∑Psd_2_ξ2_β1_Min,param.∑Psd_2_ξ2_β1_Max), (param.∑Psd_2_ξ2_β2_Min,param.∑Psd_2_ξ2_β2_Max)], NumDimensions=3, TraceMode=:silent )

			ξ1_All =  BlackBoxOptim.best_candidate(Optimization)[1]
			∑Psd_2_ξ2_β1 = BlackBoxOptim.best_candidate(Optimization)[2]
			∑Psd_2_ξ2_β2 =  BlackBoxOptim.best_candidate(Optimization)[3]
		end # if option.SubclayOpt 

		Of_AllSoil = BlackBoxOptim.best_fitness(Optimization)

		println("OPTIMAL PSD HYDRAULIC PARAMETERS:")
		println("ξ1=$ξ1_All")
		println("∑Psd_2_ξ2_β1=$∑Psd_2_ξ2_β1")
		println("∑Psd_2_ξ2_β2=$∑Psd_2_ξ2_β2")
		# println("NSE_ξ2 = $Nse_ξ2")          hay que introducir la variable en la funcion  !!!!!!!!
		println("Subclay=$Wsubclay_All")

		# Optimal parameters
		local Npsd = 100 # Just a guess
		ξ1 = zeros(Float64, Nsample)
		ξ2 = zeros(Float64, Nsample)
		Subclay = zeros(Float64, Nsample)
		Ψ_Rpart = zeros(Float64, Nsample,Npsd)
		θ_Rpart = zeros(Float64, Nsample,Npsd)
		Of_Psd = zeros(Float64, Nsample,Npsd)

		@simd for iSoil = 1:Nsample
			Subclay[iSoil] = Wsubclay_All

			ξ1[iSoil] =  ξ1_All 

			ξ2[iSoil] =  psdFunc.∑PSD_2_ξ2(∑Psd[iSoil, param.∑Psd_2_ξ2_Size]; ∑Psd_2_ξ2_β1=∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=∑Psd_2_ξ2_β2) #####################for new table model 4######
			# ξ2[iSoil] = ∑Psd_2_ξ2_β1  #####################for new table model 1######

			θ_Rpart[iSoil,1:Nrpart[iSoil]], Ψ_Rpart[iSoil,1:Nrpart[iSoil]] = psdFunc.PSD_MODEL(Rpart[iSoil,1:Nrpart[iSoil]], Psd[iSoil,:], ∑Psd[iSoil,:], Nrpart[iSoil], θsMac[iSoil], θr_Psd[iSoil], Subclay[iSoil], ξ1[iSoil], ξ2[iSoil])

			Of_Psd[iSoil] = OF_SINGLE_SOIL(ξ1[iSoil], ξ2[iSoil], Subclay[iSoil], Rpart[iSoil,1:Nrpart[iSoil]], Nrpart[iSoil], Psd[iSoil,1:Nrpart[iSoil]], ∑Psd[iSoil,1:Nrpart[iSoil]], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil], θr_Psd[iSoil])
		end # For

		return ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd
	end # function OPTIMIZATION_ALL_SOIL

#= =================================================================================================================
	       Parameters ALL SOILS
	================================================================================================================= =#
	function PARAMETERS_ALL_SOIL(Nsample, Psd, ∑Psd, Rpart, Nrpart, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)

		function OF_ALL_SOIL(ξ1, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, Nsample, Psd, ∑Psd, Rpart, Nrpart, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac, θr_Psd)
			Of_AllSoil = 0.0
			@simd for iSoil = 1:Nsample
				ξ2 =  psdFunc.∑PSD_2_ξ2(∑Psd[iSoil, param.∑Psd_2_ξ2_Size]; ∑Psd_2_ξ2_β1=param.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=param.∑Psd_2_ξ2_β2) 				

				Of_AllSoil += OF_SINGLE_SOIL(ξ1, ξ2, Subclay, Rpart[iSoil,1:Nrpart[iSoil]], Nrpart[iSoil], Psd[iSoil,1:Nrpart[iSoil]], ∑Psd[iSoil,1:Nrpart[iSoil]], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil], θr_Psd[iSoil])
			end # for loop
			return Of_AllSoil
		end # function OF_ALL_SOILS


		if option.SubclayOpt 
			
			ξ1_All =  param.ξ1
			∑Psd_2_ξ2_β1 = param.∑Psd_2_ξ2_β1
			∑Psd_2_ξ2_β2 = param.∑Psd_2_ξ2_β2
			Wsubclay_All =  param.Subclay

		elseif !(option.SubclayOpt)
			Wsubclay_All = param.Subclay

			ξ1_All =  param.ξ1
			∑Psd_2_ξ2_β1 = param.∑Psd_2_ξ2_β1
			∑Psd_2_ξ2_β2 = param.∑Psd_2_ξ2_β2
		end # if option.SubclayOpt 

		println("PSD HYDRAULIC PARAMETERS from PARAM:")
		println("ξ1=$ξ1_All")
		println("∑Psd_2_ξ2_β1=$∑Psd_2_ξ2_β1")
		println("∑Psd_2_ξ2_β2=$∑Psd_2_ξ2_β2")
		println("Subclay=$Wsubclay_All")

		# Optimal parameters
		local Npsd = 100 # Just a guess
		ξ1 = zeros(Float64, Nsample)
		ξ2 = zeros(Float64, Nsample)
		Subclay = zeros(Float64, Nsample)
		Ψ_Rpart = zeros(Float64, Nsample,Npsd)
		θ_Rpart = zeros(Float64, Nsample,Npsd)
		Of_Psd = zeros(Float64, Nsample,Npsd)

		@simd for iSoil = 1:Nsample
			Subclay[iSoil] = Wsubclay_All

			ξ1[iSoil] =  ξ1_All 

			ξ2[iSoil] =  psdFunc.∑PSD_2_ξ2(∑Psd[iSoil, param.∑Psd_2_ξ2_Size]; ∑Psd_2_ξ2_β1=param.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=param.∑Psd_2_ξ2_β2) 			

			θ_Rpart[iSoil,1:Nrpart[iSoil]], Ψ_Rpart[iSoil,1:Nrpart[iSoil]] = psdFunc.PSD_MODEL(Rpart[iSoil,1:Nrpart[iSoil]], Psd[iSoil,:], ∑Psd[iSoil,:], Nrpart[iSoil], θsMac[iSoil], θr_Psd[iSoil], Subclay[iSoil], ξ1[iSoil], ξ2[iSoil])

			Of_Psd[iSoil] = OF_SINGLE_SOIL(ξ1[iSoil], ξ2[iSoil], Subclay[iSoil], Rpart[iSoil,1:Nrpart[iSoil]], Nrpart[iSoil], Psd[iSoil,1:Nrpart[iSoil]], ∑Psd[iSoil,1:Nrpart[iSoil]], θsMac[iSoil], θr[iSoil], ΨkgMat[iSoil], σMat[iSoil], θsMat[iSoil], ΨkgMac[iSoil], σMac[iSoil], θr_Psd[iSoil])
		end # For

		return ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, ξ1, ξ2, Ψ_Rpart, θ_Rpart, Psd, ∑Psd, Of_Psd
	end # function PARAMETERS_ALL_SOIL

end # module PSD