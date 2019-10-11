# =============================================================
#		MODULE: hydroRealation
# =============================================================
module hydroRelation
	using ..param
	using BlackBoxOptim
	export OPTIMIZATION_σ_2_Ψm, Ψm_2_σ, σ_2_Ψm, OPTIMIZATION_σ_2_Ψm

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψm_2_σ(iSoil, hydro)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψm_2_σ(iSoil::Int, hydro)
			hydro.σ[iSoil] = param.hydro.Pσ_1 * ( log(hydro.Ψm[iSoil]) -1.0 ) ^ param.hydro.Pσ_2
			hydro.σ[iSoil] = max(min(hydro.σ[iSoil], param.hydro.σ_Max), param.hydro.σ_Min)
			return hydro
		end  # function: Ψm_2_σ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψm_2_σ(iSoil, hydro)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function σ_2_Ψm(iSoil::Int, hydro)
			Hydro.Ψm[iSoil] = exp(exp( inv(param.hydro.Pσ_2) * log(hydro.σ[iSoil] / param.hydro.Pσ_1)) + 1.0) #[mm]
			Hydro.Ψm[iSoil] = max(min(hydro.Ψm[iSoil] , param.hydro.Ψm_Max), param.hydro.Ψm_Min)
			return hydro
		end # function: σ_2_Ψm


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPTIMIZATION_σ_2_Ψm
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OPTIMIZATION_σ_2_Ψm(Ψm_Obs, σ_Obs)

			function OF_σ_2_Ψm(SampleTrue, Ψm_Obs, σ_Obs, Pσ_1, Pσ_2)
				Of = 0.0
				i = 1
				for σ in σ_Obs
					Ψm_Sim = exp( exp(inv(Pσ_2) * log(σ / Pσ_1)) + 1.0 ) #[mm]
					Ψm_Sim = max(min(Ψm_Sim , param.hydro.Ψm_Max), param.hydro.Ψm_Min)
					Of += (log(Ψm_Sim) - log(Ψm_Obs[i]))^2.  
					i += 1
				end # for σ in σ_Obs
				return Of 
			end # function: OF_σ_2_Ψm

			Optimization = BlackBoxOptim.bboptimize(Pσ ->  OF_σ_2_Ψm(Ψm_Obs, σ_Obs, Pσ[1], Pσ[2]); SearchRange = [(0., 100.), (0., 100.)], NumDimensions=2, TraceMode=:silent)
				
			Pσ_1 = BlackBoxOptim.best_candidate(Optimization)[1]
			Pσ_2 = BlackBoxOptim.best_candidate(Optimization)[2]

			println("Pσ_1 = $Pσ_1, Pσ_2 = $Pσ_2")
			return Pσ_1, Pσ_2
		end # function OPTIMIZATION_σ_2_Ψm
	
end  # module: hydroRealation
# ............................................................