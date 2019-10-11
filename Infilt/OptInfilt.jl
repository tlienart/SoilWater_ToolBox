
# =============================================================
#		MODULE: optInfilt
# =============================================================
module optInfilt

	# =============================================================
	#		MODULE: kg
	# =============================================================
	module kg
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPT_INFILTRATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OPT_INFILTRATION_BEST()

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
			
		end  # function: OPT_INFILTRATION
		
	end  # module: kg
	# ............................................................



end  # module: optInfilt
# ............................................................
			