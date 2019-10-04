module best
	using ..cst, ..kunsat, ..param

	using BlackBoxOptim
	export  B_2_C,  A, B, TIME_STEADYTRANSIT, B_2_Krθini,array, Q1D_TRANSIT_η,  Q1D_TRANSIT, BESTG, TIMEsteadTransit_2_SORPTIVITY


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BESTG
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BESTG(Time, θ_Ini, θs, θr, Ks, RingRadius, Sorptivity, Kr_θini; Flag_Best="Best_Gi")
			function INFILTRATION_3D_TRANSIT(Time, Sorptivity, A, B, Ks)
				return Inf = Sorptivity * (Time ^ 0.5) + (A * (Sorptivity ^ 2.0) + B * Ks) * Time
			end

			function INFILTRATION_3D_STEADY_BESTG(Time, Sorptivity, A, B, Ks)
				return Inf = (A * Sorptivity ^ 2.0 + Ks) * Time + best.C(B) * (Sorptivity^2.) / Ks
			end

			function  INFILTRATION_3D_STEADY_BESTGI(Time, Time_TransStead, Sorptivity, A, B, Ks)
				return Inf = (A * Sorptivity ^ 2.0 + Ks) * (Time - Time_TransStead) + INFILTRATION_3D_TRANSIT(Time_TransStead, Sorptivity, A, B, Ks)
			end

			B = best.B(Kr_θini)

			A = best.A(θs, θ_Ini, RingRadius)

			Time_TransStead = TIME_STEADYTRANSIT(Sorptivity, B, Ks)

			if Time <= Time_TransStead 
				return Inf = INFILTRATION_3D_TRANSIT(Time, Sorptivity, A, B, Ks)
			else
				if Flag_Best == "Best_G"
					return Inf = INFILTRATION_3D_STEADY_BESTG(Time, Sorptivity, A, B, Ks)
				else
					return Inf = INFILTRATION_3D_STEADY_BESTGI(Time, Time_TransStead, Sorptivity, A, B, Ks)
				end
			end # if
		end  # function: BESTG



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : A
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function A(θs, θ_Ini, RingRadius)
			return A = cst.γ / ( RingRadius * (θs - θ_Ini) ) # Units (mm-1)
		end # function: A


	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : B
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function B(Kr_θini)
			return B = (2.0 - cst.β) / 3.0 + Kr_θini * (1.0 + cst.β) / 3.0
		end # function: B



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : C
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function C(B)
			return C = log(1.0 / cst.β) * (1.0 + cst.β) / (6.0 * (1.0 - cst.β) * (1.0 - B) )
		end # function: C



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : B_2_C
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function B_2_C(Kr_θini)
			return C = log(1. / cst.β) / (2.0 * (1.0 -cst.β) * (1.0 - Kr_θini))
		end # function: B_2_C



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : B_2_Krθini
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function B_2_Krθini(B)
			return Kr_θini = (3.0 *B - 2.0 + cst.β) / (1.0 + cst.β)
		end # function : B_2_Krθini



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_STEADYTRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_STEADYTRANSIT(Sorptivity, B, Ks)
			return Time_TransStead = ( Sorptivity / (Ks * 2.0 * (1.0 - B)) ) ^ 2.0
		end # function: TIME_STEADYTRANSIT
	


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_1D_TRANSIT(Time, Sorptivity, B, Ks)
			return Inf = Sorptivity * (Time ^ 0.5) + B * Ks * Time
		end # function: INFILTRATION_1D_TRANSIT



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function  INFILTRATION_1D_STEADY(Time, Sorptivity, A, B, Ks)
			return Inf = Ks * Time + best.C(B) * (Sorptivity ^ 2.0) / Ks
		end # function: INFILTRATION_1D_STEADY



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Q1D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Q1D_TRANSIT(Time, Sorptivity, Se_Ini, θs, θr, σ, Ks)
			Kr_θini = kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, 1., θs, σ)
			return Q1D = 0.5 * Sorptivity * (Time^-0.5) + B(Kr_θini)*Ks
		end # Q1D_TRANSIT



    #= =============== KOSUGI =============== =#
	module kg
		include("Cst.jl")
		include("Kunsat.jl")
		include("Param.jl")
		include("Wrc.jl")
		include("Sorptivity.jl")
		include("Stat.jl")
		# include("ThetaTime.jl")
		include("HydroRelationship.jl")
		# include("QuasiExact.jl")

		using ..best
		using BlackBoxOptim, Optim
		export  BESTG_INVERSE, BESTG_INVERSE_SORPTIVITY,  STEADY_TRANSIT_BESTG,  BESTG_INVERSE_σmodel


		function STEADY_TRANSIT_BESTG(θ_Ini, θs, θr, Hkg, σ, Ks, RingRadius)
			Sorptivity = sorptivity.kg.SORPTIVITY(θ_Ini, θs, θs, θr, Hkg, σ, Ks)

			Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)

			Kr_θini= kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, 1., θs, σ)
			
			B = best.B(Kr_θini)
	
			A = best.A(θs, θ_Ini, RingRadius)
	
			return Time_TransStead_Ubest = best.TIME_STEADYTRANSIT(Sorptivity, B, Ks)
		end


		function BESTG_INVERSE(Time, iT_N, Inf_3D_Obs, θ_Ini, Se_Ini, θs, RingRadius, iT_TransStead, Time_TransStead; Option_Opt_σ=true, Flag_Best="Best_Gi", σ=1.)
			Inf_Inv = Array{Float64}(undef, iT_N)
			θr = 0.
			Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)

			P_Stead = 2.
			P_Trans = 2.
			W = 0.2
		  
			function  OF_BEST(Time, Inf_3D_Obs, iT_N, θ_Ini, θs, θr, Hkg, σ, Ks, RingRadius, iT_TransStead, Time_TransStead, Flag_Best, Se_Ini, P_Stead, P_Trans)
				Inf_Sim = zeros(iT_N)

				Sorptivity = sorptivity.kg.SORPTIVITY(θ_Ini, θs, θs, θr, Hkg, σ, Ks)

				Kr_θini= kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, 1., θs, σ)
				
				OF = 0.
				for iT in 1:iT_N
					Inf_Sim[iT]= best.BESTG(Time[iT], θ_Ini, θs, θr, Ks, RingRadius, Sorptivity, Kr_θini; Flag_Best=Flag_Best)
				end

				Of_Trans = stats.NASH_SUTCLIFFE_OF((Inf_3D_Obs[2:iT_TransStead-1]), (Inf_Sim[2:iT_TransStead-1]), P_Trans)
	
				Of_Stead =stats.NASH_SUTCLIFFE_OF(log10.(Inf_3D_Obs[iT_TransStead:iT_N]), log10.(Inf_Sim[iT_TransStead:iT_N]), P_Stead)

				Wof = W * Of_Trans + (1.0 - W) * Of_Stead
				return Wof
			end

			if Option_Opt_σ
				Optimization = BlackBoxOptim.bboptimize(X -> OF_BEST(Time[1:iT_N], Inf_3D_Obs[1:iT_N], iT_N, θ_Ini, θs, θr, 10.0 ^X[1], X[2], 10.0 ^X[3], RingRadius,iT_TransStead, Time_TransStead, Flag_Best, Se_Ini, P_Stead, P_Trans); SearchRange = [(log10(param.Hkg_Min), log10(param.Hkg_Max)), (param.σ_Min, param.σ_Max), (param.Ks_Log_Min, param.Ks_Log_Max)], NumDimensions=3, TraceMode=:silent)
				
				Hkg_Inv =10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[1])
				σ_Inv = BlackBoxOptim.best_candidate(Optimization)[2]
				Ks_Inv = 10.0 ^(BlackBoxOptim.best_candidate(Optimization)[3])
			
			else
				Optimization = BlackBoxOptim.bboptimize(X -> OF_BEST(Time[1:iT_N], Inf_3D_Obs[1:iT_N], iT_N, θ_Ini, θs, θr, 10.0 ^X[1], σ, 10.0 ^X[2], RingRadius, iT_TransStead, Time_TransStead, Flag_Best, P_Stead, P_Trans); SearchRange = [(log10(param.Hkg_Min), log10(param.Hkg_Max)), (param.Ks_Log_Min, param.Ks_Log_Max)], NumDimensions=2, TraceMode=:silent)
				Hkg_Inv =10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[1])
				Ks_Inv = 10.0 ^ (BlackBoxOptim.best_candidate(Optimization)[2])
				σ_Inv = σ
			end
			
			Kr_θini_Inv = kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ_Inv, 1., θs, σ_Inv)

			Sorptivity_Inv = sorptivity.kg.SORPTIVITY(θ_Ini, θs, θs, θr, Hkg_Inv, σ_Inv, Ks_Inv)

			# Computing best Inf_Inv
			for iT in 1: iT_N
				Inf_Inv[iT] = best.BESTG(Time[iT], θ_Ini, θs, θr, Ks_Inv, RingRadius, Sorptivity_Inv, Kr_θini_Inv)
			end

			return Sorptivity_Inv, Hkg_Inv, σ_Inv, Ks_Inv, Kr_θini_Inv, Inf_Inv
		end



		function BESTG_INVERSE_σmodel(Sorpt, Ks, Time, iT_N, Inf_3D_Obs, θ_Ini, Se_Ini, θs, RingRadius, iT_TransStead, Time_TransStead; Flag_Best="Best_Gi")
			Inf_Inv = Array{Float64}(undef, iT_N)
			θr = 0.
			Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)
			P_Trans = 2.
			P_Stead = 2.
			W = 0.2
		  
			function  OF_BEST_σmodel(Sorpt, Time, Inf_3D_Obs, iT_N, θ_Ini, θs, θr, σ, Ks, RingRadius, iT_TransStead, Time_TransStead, Flag_Best, Se_Ini, P_Stead, P_Trans)
				Inf_Sim = zeros(iT_N)

				Kr_θini= kunsat.kg.KUNSAT(Se_Ini, θs, θr, σ, 1., θs, σ)

					
				for iT in 1:iT_N
					Inf_Sim[iT]= best.BESTG(Time[iT], θ_Ini, θs, θr, Ks, RingRadius, Sorpt, Kr_θini; Flag_Best=Flag_Best)
				end

				Of_Trans = stats.NASH_SUTCLIFFE_OF((Inf_3D_Obs[2:iT_TransStead-1]), (Inf_Sim[2:iT_TransStead-1]), P_Trans)
	
				Of_Stead =stats.NASH_SUTCLIFFE_OF(log10.(Inf_3D_Obs[iT_TransStead:iT_N]), log10.(Inf_Sim[iT_TransStead:iT_N]), P_Stead)
		
				Hkg_Sorpt = sorptivity.kg.SORPTIVITY_2_Hkg(Sorpt, θ_Ini, θs, θr, σ, Ks)

				Hkg_σ = relationship.σ_2_Hkg(σ)
				
				Of_Hkg = abs(log10(Hkg_Sorpt) - log10(Hkg_σ)) / 10.
				
				Wof = W * Of_Trans + (1.0 - W) * Of_Stead + Of_Hkg 

				return Wof
			end

			# Optimization = BlackBoxOptim.bboptimize(X -> OF_BEST_σmodel(Sorpt, Time[1:iT_N], Inf_3D_Obs[1:iT_N], iT_N, θ_Ini, θs, θr, X[1], Ks, RingRadius, iT_TransStead, Time_TransStead, Flag_Best, Se_Ini, P_Stead, P_Trans); SearchRange = [(param.σ_Min, param.σ_Max)], NumDimensions=1, TraceMode=:silent)

			Optimization =  Optim.optimize(σ -> OF_BEST_σmodel(Sorpt, Time[1:iT_N], Inf_3D_Obs[1:iT_N], iT_N, θ_Ini, θs, θr, σ, Ks, RingRadius, iT_TransStead, Time_TransStead, Flag_Best, Se_Ini, P_Stead, P_Trans), param.σ_Min, param.σ_Max, GoldenSection() )
            
			σ = Optim.minimizer(Optimization)[1]

			Hkg_Sorpt = sorptivity.kg.SORPTIVITY_2_Hkg(Sorpt, θ_Ini, θs, θr, σ, Ks)
			
			return σ, Hkg_Sorpt
		end

		

		function BESTG_INVERSE_SORPTIVITY(Se_Ini_2, Time, Inf_3D_Obs, Se_Ini, θ_Ini, θs, θr, Ks, RingRadius, Time_TransStead, iT_TransStead, iT_N, Sorptivity_1)
			Se_Ini_3 = 0.7

			θ_Ini_2 = wrc.se.Se_2_θ(Se_Ini_2, θs, θr)
			θ_Ini_3 = wrc.se.Se_2_θ(Se_Ini_3, θs, θr)
			θr = 0.
			
			function OF_BEST(Time, Inf_3D_Obs, RingRadius, θ_Ini, θs, θr, σ, Ks, iT_TransStead, iT_N, Sorptivity_1, Se_Ini_2, θ_Ini_2)
				Inf_Sim = Array{Float64}(undef, iT_N)

				# For θ_Ini_2
				
				θ_Time, Of_θ = thetaTime.kg.θTIME(Time[1:iT_TransStead], Inf_3D_Obs[1:iT_TransStead], RingRadius, Sorptivity_1, θ_Ini, θs, θr, σ, Ks, iT_TransStead, iT_N)

				iθ_Ini_2 = array.SEARCH_INDEX(θ_Time[1:iT_TransStead], θ_Ini_2) # getting the closest index
				iθ_Ini_3 = array.SEARCH_INDEX(θ_Time[1:iT_TransStead], θ_Ini_2) # getting the closest index

				Hkg = sorptivity.kg.SORPTIVITY_2_Hkg(Sorptivity_1, θ_Ini, θs, θr, σ, Ks)

				Kr_θini_2 = kunsat.kg.KUNSAT(Se_Ini_2, θs, θr, σ, 1., θs, σ)
				Kr_θini_3 = kunsat.kg.KUNSAT(Se_Ini_3, θs, θr, σ, 1., θs, σ)
				
				Sorptivity_2 = sorptivity.kg.SORPTIVITY(θ_Ini_2, θs, θs, θr, Hkg, σ, Ks)
				Sorptivity_3 = sorptivity.kg.SORPTIVITY(θ_Ini_3, θs, θs, θr, Hkg, σ, Ks)
				

				# For θ_Ini_2
				Inf_Sim[iθ_Ini_2-1] =  Inf_3D_Obs[iθ_Ini_2-1]
				Wof = 0.	
				for iT in iθ_Ini_2:iT_N
					Inf_Sim[iT] = best.BESTG(Time[iT]-Time[iθ_Ini_2-1], θ_Ini_2, θs, θr, Ks, RingRadius, Sorptivity_2, Kr_θini_2; Flag_ModSTeadyTrans=true) +  Inf_3D_Obs[iθ_Ini_2-1]
				end

				Of_Trans = stats.NASH_SUTCLIFFE_OF(Inf_3D_Obs[iθ_Ini_2:iT_TransStead], Inf_Sim[iθ_Ini_2:iT_TransStead], 2.)

				Of_Stead =stats.NASH_SUTCLIFFE_OF(Inf_3D_Obs[iT_TransStead-1:iT_N], Inf_Sim[iT_TransStead-1:iT_N], 1.)

				Wof =  Wof + 0.5 * Of_Trans + 0.5 * Of_Stead + Of_θ

				# println("σ= $σ,  Hkg=$Hkg,  Of_Trans=$Of_Trans,  Of_Stead=$Of_Stead" )

				return Wof, Inf_Sim
			end #OF_BEST	


			Optimization = BlackBoxOptim.bboptimize(σ  -> OF_BEST(Time[1:iT_N], Inf_3D_Obs[1:iT_N], RingRadius, θ_Ini, θs, θr, σ[1], Ks, iT_TransStead, iT_N, Sorptivity_1, Se_Ini_2, θ_Ini_2)[1]; SearchRange = [(param.σ_Min, param.σ_Max)], NumDimensions=1, TraceMode=:silent)
			
			σ = BlackBoxOptim.best_candidate(Optimization)[1]
		
			~, Inf_Sim = OF_BEST(Time[1:iT_N], Inf_3D_Obs[1:iT_N], RingRadius, θ_Ini, θs, θr, σ, Ks, iT_TransStead, iT_N, Sorptivity_1, Se_Ini_2, θ_Ini_2)

			Hkg = sorptivity.kg.SORPTIVITY_2_Hkg(Sorptivity_1, θ_Ini, θs, θr, σ, Ks)

			println( "Se_Ini_2=  ", Se_Ini_2, "  ,σ=" , σ, "  ,Hkg=", Hkg,  "\n ")

			return Inf_Sim, σ, Hkg, Ks
		end #BESTG_INVERSE_SORPTIVITY



	end # Module Kosugi
	
    #= =============== VAN GENUCHTEN =============== =#
	module vg
			# include("Cst.jl")
			include("Kunsat.jl")
			include("Param.jl")
			include("Wrc.jl")
			include("Sorptivity.jl")
			include("Stat.jl")
			# include("ThetaTime.jl")
			# include("HydroRelationship.jl")
			# include("QuasiExact.jl")
	using BlackBoxOptim
	using ..best
	export BEST_INVERSE_SORPTIVITY, BEST_INVERSE_IMPROVED, STEADY_TRANSIT_BESTG

		function STEADY_TRANSIT_BESTG(θ_Ini, θs, θr, Hvg, N, Ks, Km, RingRadius)
			Sorptivity = sorptivity.vg.SORPTIVITY(θ_Ini, θs, θs, θr, Hvg, N, Ks, Km)
		
			Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)

			Kr_θini= kunsat.vg.KUNSAT(Se_Ini, N, 1., Km)
			
			B = best.B(Kr_θini)

			A = best.A(θs, θ_Ini, RingRadius)

			return Time_TransStead_Ubest = best.TIME_STEADYTRANSIT(Sorptivity, B, Ks)
		end



		function BESTG_INVERSE(Time, iT_N, Inf_3D_Obs, θ_Ini, Se_Ini, θs, Km, RingRadius, iT_TransStead, Time_TransStead; Option_Opt_N=true, Flag_Best="Best_Gi", N=1.5)
			Inf_Inv = Array{Float64}(undef, iT_N)
			
			θr = 0.

			function  OF_BEST(Time, Inf_3D_Obs, iT_N, θ_Ini, θs, θr, Hvg, N, Ks, RingRadius, Km)
				Inf_Sim = zeros(iT_N)

				Sorptivity = sorptivity.vg.SORPTIVITY(θ_Ini, θs, θs, θr, Hvg, N, Ks, Km)

				Se_Ini = wrc.se.θ_2_Se(θ_Ini, θs, θr)

				Kr_θini= kunsat.vg.KUNSAT(Se_Ini, N, 1., Km)
				
				OF = 0.
				for iT in 1:iT_N
					Inf_Sim[iT]= best.BESTG(Time[iT], θ_Ini, θs, θr, Ks, RingRadius, Sorptivity, Kr_θini; Flag_Best=Flag_Best)
				end

				OF_Trans = stats.NASH_SUTCLIFFE_OF(Inf_3D_Obs[2:iT_TransStead], (Inf_Sim[2:iT_TransStead]), 4.)
		
				OF_Stead =stats.NASH_SUTCLIFFE_OF((Inf_3D_Obs[iT_TransStead:iT_N]), (Inf_Sim[iT_TransStead:iT_N]), 1.)
		
				OF = OF_Trans + OF_Stead
				return OF
			end

			if Option_Opt_N
				# N_Min and N_Max depends on Km
				if Km == 1
					N_Min = param.N_Km1_Min
					N_Max = param.N_Km1_Max
				elseif  Km == 2
					N_Min = param.N_Km2_Min
					N_Max = param.N_Km2_Max
				end

				Optimization = bboptimize(X -> OF_BEST(Time[1:iT_N], Inf_3D_Obs[1:iT_N], iT_N, θ_Ini, θs, θr, 10.0 ^X[1], X[2], 10.0 ^X[3], RingRadius, Km); SearchRange = [(log10(param.Hvg_Min), log10(param.Hvg_Max)), (N_Min, N_Max), (param.Ks_Log_Min, param.Ks_Log_Max)], NumDimensions=3, TraceMode=:silent)

					Hvg_Inv =10.0 ^(best_candidate(Optimization)[1])
					N_Inv = best_candidate(Optimization)[2]
					Ks_Inv = 10. ^(best_candidate(Optimization)[3])
			
				else
				Optimization = bboptimize(X -> OF_BEST(Time[1:iT_N], Inf_3D_Obs[1:iT_N], iT_N, θ_Ini, θs, θr, 10.0 ^X[1], N, 10.0 ^X[2], RingRadius, Km); SearchRange = [(log10(param.Hvg_Min), log10(param.Hvg_Max)), (param.Ks_Log_Min, param.Ks_Log_Max)], NumDimensions=2, TraceMode=:silent)
					Hvg_Inv =10.0 ^(best_candidate(Optimization)[1])
					Ks_Inv = 10. ^(best_candidate(Optimization)[2])
					N = N_Inv
			end
			
			Kr_θini_Inv= kunsat.vg.KUNSAT(Se_Ini, N_Inv, 1., Km)
			
			Sorptivity_Inv = sorptivity.vg.SORPTIVITY(θ_Ini, θs, θs, θr, Hvg_Inv, N_Inv, Ks_Inv, Km)

			# Computing best Inf_Inv
			for iT in 1: iT_N
				Inf_Inv[iT] = best.BESTG(Time[iT], θ_Ini, θs, θr, Ks_Inv, RingRadius, Sorptivity_Inv, Kr_θini_Inv)
			end

			return Sorptivity_Inv, Hvg_Inv, N_Inv, Ks_Inv, Kr_θini_Inv, Inf_Inv
		end

	end # module vg

end # module best
