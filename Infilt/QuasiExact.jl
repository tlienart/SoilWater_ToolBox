module quasiExact # quasi-exact objective function
	import ..option, ..sorptivity, ..wrc, ..kunsat, ..option, ..param
	import BlackBoxOptim, Optim
   	export INFILTRATION3D_2_1D, OF_INFILTRATION_2_HYDRO


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_2_TIMEη
	#		= COMPUTE NORMALISED TIME: Tη =#
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_2_TIMEη(iT, Sorptivity, ΔK)
			return Time_η = iT * 2.0 * (ΔK / Sorptivity) ^ 2.0
		end # function: TIME_2_TIMEη


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION3D_2_1D
	# 		TRANSFORMS INFILTRATION_3D TO INFILTRATION_1D
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION3D_2_1D(iT, Inf_3D_Obs, Sorptivity, infiltParam, θs, θr)
			Δθ = θs - θr
			return Inf_1D = max(Inf_3D_Obs - (iT * cst.γ * Sorptivity ^ 2.0) / (infiltParam.RingRadius[iSoil] * Δθ), cst.ϵ)
		end # function : INFILTRATION3D_2_1D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTTRATIONη_2_INFILTRATION3D
	# 		TRANSFORMS NORMALIZED INFILTRATION TO INFILTRATION-3D
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTTRATIONη_2_INFILTRATION3D(iT, Time_η, Inf_η,Sorptivity, ΔK,  K_θini, Δθ, infiltParam)
			# Function compute infiltration-1d from normalized infiltration
			function INFILTTRATIONη_2_INFILTRATION1D(iT, Inf_η,Sorptivity, ΔK,  K_θini)
				return INFILTRATION_1D = K_θini*iT + Inf_η * (Sorptivity^2.) / (2.0 *ΔK)
			end
			ΔI_η = Time_η * cst.γ
			return Inf_3D_Obs = INFILTTRATIONη_2_INFILTRATION1D(iT, Inf_η,Sorptivity, ΔK, K_θini) + ΔI_η * (Sorptivity ^ 4.0) / (2.0 * infiltParam.RingRadius[iSoil] * Δθ *(ΔK ^ 2.0))
		end # function : INFILTTRATIONη_2_INFILTRATION3D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION3D_2_INFILTTRATIONη
   	# 		TRANSFORMS INFILTRATION-3D TO NORMALIZED INFILTRATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION3D_2_INFILTTRATIONη(iT,Inf_3D_Obs,Sorptivity, ΔK, K_θini, Δθ, infiltParam)
			Inf_η = max((2.0 *ΔK /Sorptivity^2.) * (Inf_3D_Obs  - K_θini*iT -  cst.γ * TIME_2_TIMEη(iT,Sorptivity, ΔK) * (Sorptivity^4.) /  (infiltParam.RingRadius[iSoil] * Δθ * 2.0* (ΔK^2.)) ), cst.ϵ)
			return Inf_η
		end


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_QUASIEXACTη
   	# 		NORMALISED QUASI-EXACT 3D CUMULATIVE INFILTRATION EQUATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_QUASIEXACTη(Time_η, Inf_η)
			Left_Term = Time_η
			Right_Term =  (1.0 / (1.0 -cst.β)) * (Inf_η - log((exp(cst.β*Inf_η) + cst.β-1.) / cst.β))
			if  Right_Term  < 0.
				# OF = 100000.0 *(abs(exp(cst.β*Inf_η) + cst.β-1.))
				OF = 10000.0 * exp(Inf_η)
			else
				OF =  abs(Left_Term - Right_Term)
			end
			return OF
		end # function OF_QUASIEXACTη


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_INFILTRATION_2_HYDRO
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_INFILTRATION_2_HYDRO(N_Infilt, iT_TransStead, Time, Inf_3D_Obs,Sorptivity, ΔK, K_θini, Δθ, infiltParam)
			Left_Term = zeros(N_Infilt[iSoil])
			Right_Term = zeros(N_Infilt[iSoil])

			Of_Penalty = 0.
			Of_Stead = 0.
			Of_Trans = 0.
			W = 0.2
			for iT in 2:N_Infilt[iSoil]
				Time_η = quasiExact.TIME_2_TIMEη(Time[iT],Sorptivity, ΔK)
				Inf_η = quasiExact.INFILTRATION3D_2_INFILTTRATIONη(Time[iT], Inf_3D_Obs[iT],Sorptivity, ΔK, K_θini, Δθ, infiltParam)

				Left_Term[iT] = Time_η
				Right_Term[iT] =  (1.0 / (1.0 -cst.β)) * (Inf_η - log((exp(cst.β*Inf_η) + cst.β-1.) / cst.β))
				
				if Right_Term[iT] > 0.
					if iT <= iT_TransStead
						Of_Trans = Of_Trans + ((Left_Term[iT]) - (Right_Term[iT]))^2.
					else
						Of_Stead = Of_Stead + (log10(Left_Term[iT]) - log10(Right_Term[iT]))^2.
					end
				else
					Of_Penalty = Of_Penalty + 1000.0* exp(Inf_η)
					Right_Term[iT] = 0.
				end
			end
			Of = W * Of_Trans / Float64(iT_TransStead-1)  + (1.0 -W) * Of_Stead / Float64(N_Infilt[iSoil] - iT_TransStead + 1) + Of_Penalty

			return Of
		end # function:  OF_INFILTRATION_2_HYDRO


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDRO_2_INFILTRATION3D
	# 		COMPUTE INFILTRATION_3D FROM OPTIMIZED HYDRAULIC PARAMETERS
    # 		Solving quasiexact solution
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
		function HYDRO_2_INFILTRATION3D(Time,Sorptivity, Ks, K_θini, θ_Ini, N_Infilt, iT_TransStead, hydro, iSoil, infiltParam)
			Inf_3D_Sim= Array{Float64}(undef, N_Infilt[iSoil]) # preparing the matrix
			Inf_η= Array{Float64}(undef, N_Infilt[iSoil]) # preparing the matrix
			Δθ = hydro.θs[iSoil] - θ_Ini
			ΔK = Ks - K_θini
		
			# At t=1
			Inf_3D_Sim[1] = 0.
			Inf_η[1] = 0.
			Inf_η_Min = 0.0001
			Inf_η_Max = 0.05 * Time[2] #Since Time[1] = 0

			for iT in 2:N_Infilt[iSoil] # Looping for every time step
				Time_η= TIME_2_TIMEη(Time[iT],Sorptivity, ΔK)

				# We are solving for Inf_η
				Optimization =  Optim.optimize(Inf_η ->  OF_QUASIEXACTη(Time_η, Inf_η), Inf_η_Min, Inf_η_Max, GoldenSection())
				Inf_η[iT] = Optim.minimizer(Optimization)[1]

				# Deriving the new bounds such that infiltration increases with time & the slope decreases with time
				Inf_η_Min = Inf_η[iT] + 0.0001

				# Maximum infiltration rate for T+1: (Inf[T2] - Inf[T1]) / (T2 - T1) which is 1 seconds
				if iT <= N_Infilt[iSoil] - 1
					Inf_η_Max = Inf_η[iT] + (Time[iT+1]- Time[iT]) * (Inf_η[iT] - Inf_η[iT-1]) / (Time[iT] - Time[iT-1])
				else
					Inf_η_Max = Inf_η[iT] + (Inf_η[iT] - Inf_η[iT-1])
				end

				# Transforming INFILTRATION3D 2 INFILTTRATIONη
				Inf_3D_Sim[iT] =  INFILTTRATIONη_2_INFILTRATION3D(Time[iT], Time_η, Inf_η[iT],Sorptivity, ΔK, K_θini, Δθ, infiltParam)
			end
			return Inf_3D_Sim
		end # function: HYDRO_2_INFILTRATION3D


    #= =============== KOSUGI =============== =#
	module kg	
		import ..option, ..sorptivity, ..wrc, ..kunsat, ..option, ..param
		import ..quasiExact
        export BlackBoxOptim, Optim
        export INFILTRATION3D_2_HYDRO, INFILTRATION3D_2_HYDRO_SORPTIVITY, INFILTRATION3D_2_HYDRO_σMOD

        # OPTIMIZE THE HYDRAULIC PARAMETERS FROM INFILTRATION DATA============================
        function INFILTRATION3D_2_HYDRO_σMOD(Time, Inf_3D_Obs, N_Infilt, θs, θ_Ini, iT_TransStead, Time_TransStead, Ks,Sorptivity, infiltParam)
            θr = 0.
            Δθ = θs - θ_Ini
            Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)

            # Objective function which matches observed with simulated infiltration
            function OF_Fast_η(Time, Inf_3D_Obs, Δθ, N_Infilt, infiltParam, θr, θ_Ini, σ, iT_TransStead)

                Hkg_σ = relationship.σ_2_Hkg(σ)

                Hkg_Sorpt = sorptivity.kg.SORPTIVITY_2_Hkg(Sorptivity, θ_Ini, θs, θr, σ, Ks)

                K_θini = kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, Ks, θs, σ)

                ΔK = Ks - K_θini

                Of_Hkg = abs(log10(Hkg_Sorpt) - log10(Hkg_σ))

                Wof = quasiExact.OF_INFILTRATION_2_HYDRO(N_Infilt, iT_TransStead, Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil]],Sorptivity, ΔK, K_θini, Δθ, infiltParam) + Of_Hkg / 10.

                return Wof
            end #  function INFILTRATION3D_2_HYDRO_σMOD

            # OPTIMIZATION

                Optimization =  Optim.optimize(σ ->  OF_Fast_η(Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil] ], Δθ, N_Infilt , infiltParam, θr, θ_Ini, σ, iT_TransStead), param.σ_Min, param.σ_Max, GoldenSection() )

                # Values of the optimal hydraulic params
                σ = Optim.minimizer(Optimization)[1]
            
          

            Hkg_Sorpt = sorptivity.kg.SORPTIVITY_2_Hkg(Sorptivity, θ_Ini, θs, θr, σ, Ks)
           
            return σ, Hkg_Sorpt
        end


            # OPTIMIZE THE HYDRAULIC PARAMETERS FROM INFILTRATION DATA============================
            function INFILTRATION3D_2_HYDRO(Time, Inf_3D_Obs, N_Infilt, θs, θ_Ini, infiltParam, iT_TransStead, Time_TransStead, Option_Opt_σ, σ=1.)
                θr = 0.
                Δθ = θs - θ_Ini
    
                # Objective function which matches observed with simulated infiltration
                function OF_Fast_η(Time, Inf_3D_Obs, Δθ, N_Infilt, infiltParam, θr, θ_Ini, Hkg, Ks, σ, iT_TransStead)
    
                   Sorptivity = sorptivity.kg.SORPTIVITY(θ_Ini, θs, θs, θr, Hkg, σ, Ks)
    
                    Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)
    
                    K_θini = kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, Ks, θs, σ)
    
                    ΔK = Ks - K_θini
    
                    Wof = quasiExact.OF_INFILTRATION_2_HYDRO(N_Infilt, iT_TransStead, Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil]],Sorptivity, ΔK, K_θini, Δθ,infiltParam)
    
                    Of_Sorpt = 0.
                    # Of_Sorpt = abs(log10(Hkg) - log10(relationship.σ_2_Hkg(σ)))
    
                    return Wof =  Wof +  Of_Sorpt
                end
    
                # OPTIMIZATION
                if Option_Opt_σ # If Ks is not known
                    Optimization = BlackBoxOptim.bboptimize(Param ->  OF_Fast_η(Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil] ], Δθ, N_Infilt, infiltParam, θr, θ_Ini, 10.0 ^Param[1], 10.0 ^Param[2], Param[3], iT_TransStead) ; SearchRange =[ (log10(param.Hkg_Min), log10(param.Hkg_Max)), (param.Ks_Log_Min, param.Ks_Log_Max), (param.σ_Min, param.σ_Max)], NumDimensions=3, TraceMode=:silent)
    
                    # Values of the optimal hydraulic params
                    Hkg = 10.0 ^(BlackBoxOptim.best_candidate(Optimization)[1])
                    Ks = 10.0 ^(BlackBoxOptim.best_candidate(Optimization)[2])
                    σ = BlackBoxOptim.best_candidate(Optimization)[3]
                
                else # If σ is known
                    Optimization = BlackBoxOptim.bboptimize(Param ->  OF_Fast_η(Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil]], Δθ, N_Infilt, infiltParam, θr, θ_Ini, 10.0 ^Param[1], 10.0 ^Param[2], σ, iT_TransStead) ; SearchRange = [ (log10(param.Hkg_Min), log10(param.Hkg_Max)), (param.Ks_Log_Min, param.Ks_Log_Max)], NumDimensions=2, TraceMode=:silent)
                    # Values of the optimal hydraulic params
                    Hkg = 10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[1])
                    Ks = 10.0 ^(BlackBoxOptim.best_candidate(Optimization)[2])
                end
    
               Sorptivity = sorptivity.kg.SORPTIVITY(θ_Ini, θs, θs, θr, Hkg, σ, Ks)
    
                # Converting θ_Ini to Kr_θini
                Se_Ini= wrc.se.θ_2_Se(θ_Ini, θs, θr)
                K_θini= kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, Ks, θs, σ)
                Kr_θini = K_θini / Ks
                
                return Ks, Kr_θini,Sorptivity, σ, Hkg
            end





        function INFILTRATION3D_2_HYDRO_SORPTIVITY(Time, Inf_3D_Obs, Se_Ini, θ_Ini, θs, θr, Ks, infiltParam, Time_TransStead, iT_TransStead, N_Infilt, Sorptivity_1)

            θr = 0.

            Se_Ini_2 = 0.3
            Se_Ini_3 = 0.6
            Se_Ini_4 = 0.6

            θ_Ini_2 = wrc.se.Se_2_θ(Se_Ini_2, θs, θr)
            θ_Ini_3 = wrc.se.Se_2_θ(Se_Ini_3, θs, θr)
            θ_Ini_4 = wrc.se.Se_2_θ(Se_Ini_4, θs, θr)
            
            function OF_TIME(Time, Inf_3D_Obs, infiltParam, θ_Ini, θs, θr, σ, Ks, iT_TransStead, N_Infilt, Sorptivity_1, Se_Ini_2, θ_Ini_2,  Se_Ini_3 ,θ_Ini_3)

                Of_θ, θ_Time = thetaTime.kg.θTIME(Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil]], InfiltParam, Sorptivity_1, θ_Ini, θs, θr, σ, Ks, iT_TransStead, N_Infilt[iSoil])
                Wof = Of_θ

                return Wof, θ_Time
            end #OF_BEST
            
            function OF_θ_Ini(Time, Inf_3D_Obs, infiltParam, θ_Ini, θs, θr, σ, Ks, iT_TransStead, N_Infilt, Sorptivity_1, Se_Ini_2, θ_Ini_2,  Se_Ini_3 ,θ_Ini_3)

					 Of_θ, θ_Time = thetaTime.kg.θTIME(Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil]], infiltParam, Sorptivity_1, θ_Ini, θs, θr, σ, Ks, iT_TransStead, N_Infilt)
					 

              
                iθ_Ini_4 = array.SEARCH_INDEX(θ_Time[1:N_Infilt[iSoil]], θ_Ini_4) # getting the closest index

                Hkg = sorptivity.kg.SORPTIVITY_2_Hkg(Sorptivity_1, θ_Ini, θs, θr, σ, Ks)

                K_θini_4 = kunsat.kg.KUNSAT(Se_Ini_4, θs, θr, σ, Ks, θs, σ)
                
                Sorptivity_4 = sorptivity.kg.SORPTIVITY(θ_Ini_4, θs, θs, θr, Hkg, σ, Ks)
    
                ΔK_4 = Ks - K_θini_4
                Δθ_4 = θs - θ_Ini_4
                Wof_4 = quasiExact.OF_INFILTRATION_2_HYDRO(N_Infilt[iSoil]-iθ_Ini_4+1, iT_TransStead-iθ_Ini_4+1, Time[iθ_Ini_4:N_Infilt[iSoil]]-Time[iθ_Ini_4-1], Inf_3D_Obs[iθ_Ini_4:N_Infilt[iSoil]]-Inf_3D_Obs[iθ_Ini_4-1], Sorptivity_4, ΔK_4, K_θini_4, Δθ_4, infiltParam)

					 Wof = Wof_4 / (iT_TransStead-iθ_Ini_4+1)
					#  println("$Of_θ  $Wof_4")
					 return Wof
            end
                         


            Optimization = BlackBoxOptim.bboptimize(σ  ->  OF_TIME(Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil]], infiltParam, θ_Ini, θs, θr, σ[1], Ks, iT_TransStead, N_Infilt[iSoil], Sorptivity_1, Se_Ini_2, θ_Ini_2, Se_Ini_3, θ_Ini_3)[1]; SearchRange = [(param.σ_Min, param.σ_Max)], NumDimensions=1, TraceMode=:silent)
            
            σ = BlackBoxOptim.best_candidate(Optimization)[1]
        
            Hkg = sorptivity.kg.SORPTIVITY_2_Hkg(Sorptivity_1, θ_Ini, θs, θr, σ, Ks)

            Wof, θ_Time = OF_TIME(Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil]], infiltParam, θ_Ini, θs, θr, σ, Ks, iT_TransStead, N_Infilt[iSoil], Sorptivity_1, Se_Ini_2, θ_Ini_2, Se_Ini_3, θ_Ini_3)

            println( "Se_Ini_2=  ", Se_Ini_2, "  ,σ=" , σ, "  ,Hkg=", Hkg,  "\n ")
            println("OF= $Wof")

            return σ, Hkg, Ks, θ_Time
        end #BESTG_INVERSE_SORPTIVITY

    end


    #= =============== VAN GENUCHTEN =============== =#
	module vg
	include("Cst.jl")
	include("Param.jl")
	include("Wrc.jl")
	include("Kunsat.jl")
	include("Sorptivity.jl")
	include("Stat.jl")
	include("Tools.jl")
	include("HydroRelationship.jl")
	using ..quasiExact
    using  BlackBoxOptim, Optim
        export INFILTRATION3D_2_HYDRO

        # OPTIMIZE THE HYDRAULIC PARAMETERS FROM INFILTRATION DATA============================
        function INFILTRATION3D_2_HYDRO(Time, Inf_3D_Obs, N_Infilt, θs, θ_Ini, infiltParam, Km, iT_TransStead, Time_TransStead, Option_Opt_N, N=1.)
            θr = 0.
            Δθ = θs - θ_Ini
            
            # Objective function which matches observed with simulated infiltration
            function OF_Fast_η(Time, Inf_3D_Obs, Δθ, N_Infilt, infiltParam, θr, θ_Ini, Hvg, Ks, N, Km)
            
               Sorptivity = sorptivity.vg.SORPTIVITY(θ_Ini, θs, θs,  θr, Hvg, N, Ks, Km)

                Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)
            
                K_θini = kunsat.vg.KUNSAT(Se_Ini, N, Ks, Km)
            
                ΔK = Ks - K_θini

                OF_Cumul = quasiExact.OF_INFILTRATION_2_HYDRO(N_Infilt[iSoil], iT_TransStead, Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil]],Sorptivity, ΔK, K_θini, Δθ, infiltParam)
                return OF_Cumul
            end

            # OPTIMIZATION
            # N_Min and N_Max depends on Km
            if Km == 1
                N_Min = param.N_Km1_Min
                N_Max = param.N_Km1_Max
            elseif  Km == 2
                N_Min = param.N_Km2_Min
                N_Max = param.N_Km2_Max
            end

            if Option_Opt_N # If Ks is not known
                Optimization = BlackBoxOptim.bboptimize(Param ->  OF_Fast_η(Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil]], Δθ, N_Infilt[iSoil], infiltParam,  θr, θ_Ini, 10.0 ^Param[1], 10.0 ^Param[2], Param[3], Km) ; SearchRange =[ (log10(param.Hvg_Min), log10(param.Hvg_Max)), (log10(param.Ks_Min), log10(param.Ks_Max)), (N_Min, N_Max)], NumDimensions=3, TraceMode=:silent)
                # Values of the optimal hydraulic params
                Hvg = 10.0 ^(BlackBoxOptim.best_candidate(Optimization)[1])
                Ks = 10.0 ^(BlackBoxOptim.best_candidate(Optimization)[2])
                N = BlackBoxOptim.best_candidate(Optimization)[3]
            else # if N is known
                Optimization = BlackBoxOptim.bboptimize(Param ->  OF_Fast_η(Time[1:N_Infilt[iSoil]], Inf_3D_Obs[1:N_Infilt[iSoil]], Δθ, N_Infilt[iSoil], infiltParam,  θr, θ_Ini, 10.0 ^Param[1], 10.0 ^Param[2], N, Km) ; SearchRange = [ (log10(param.Hvg_Min), log10(param.Hvg_Max)), (param.Ks_Log_Min, param.Ks_Log_Max)], NumDimensions=2, TraceMode=:silent)
                # Values of the optimal hydraulic params
                Hvg = 10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[1])
                Ks = 10.0 ^(BlackBoxOptim.best_candidate(Optimization)[2])      
            end

           Sorptivity = sorptivity.vg.SORPTIVITY(θ_Ini, θs, θs, θr, Hvg, N, Ks, Km)

            Se_Ini= wrc.se.θ_2_Se(θ_Ini, θs, θr)
 
            Kr_θini=  kunsat.vg.KUNSAT(Se_Ini, N, 1., Km)
            
            return Ks, Kr_θini,Sorptivity, N, Hvg
        end
    end

end