# =============================================================
#		module: startKsModel
# =============================================================
module startKsModel
	import ..θψ2Ks, ..stats
	export START_KSMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_KSMODEL(hydro, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], OptLayer=1)

			# Just run the model
			if optimKsmodel.Flag_Opt
				SearchRange = SEARCHRANGE(optimKsmodel)
				println(SearchRange)

				Optimization = BlackBoxOptim.bboptimize(X -> OF_KSMODEL(hydro, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel, optimKsmodel, OptLayer, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[]); SearchRange=SearchRange, NumDimensions=optimKsmodel.NparamOpt, TraceMode=:silent)

				X = BlackBoxOptim.best_candidate(Optimization)

				ksmodelτ = X_2_τ(ksmodelτ, optimKsmodel, OptLayer, X)

				Kₛ_Model = KS_MODEL(hydro, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

			else
				Kₛ_Model = KS_MODEL(hydro, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])
	
			end  # if: optimKsmodel
			
			for iZ=1:N_iZ
				if "Ks" ∉ optim.ParamOpt
					hydro.Ks[iZ] = Kₛ_Model[iZ]
				end #  hydro.Ks[iZ] < eps(100.0)
			end # if: hydro.Ks[iZ] > eps(10.0)
			
		return hydro, Kₛ_Model
		end  # function: START_KSMODEL


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KS_MODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function KS_MODEL(hydro, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

				for iZ=1:N_iZ
					if Flag_RockFragment
						RockFragment₁ = RockFragment[iZ]
					else
						RockFragment₁ = 0.0
					end #@isdefined RockFragment
					if Flag_IsTopsoil
						IsTopsoil₁ = Int64(IsTopsoil[iZ])
					else
						IsTopsoil₁ = 1					
					end  # if: @isdefined IsTopsoil

					Kₛ_Model[iZ] = θψ2Ks.θΨ_2_KS(hydro, IsTopsoil₁, iZ, ksmodelτ; RockFragment=RockFragment₁)

					@show Kₛ_Model[iZ] 
				end # if: hydro.Ks[iZ] > eps(10.0)	
			return Kₛ_Model
			end  # function: KS_MODEL


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : OF_KSMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function OF_KSMODEL(hydro, Kₛ_Model, ksmodelτ, N_iZ, optim, optimksmodel, optimKsmodel, OptLayer, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

				ksmodelτ = X_2_τ(ksmodelτ, optimksmodel, OptLayer, X)

				Kₛ_Model = KS_MODEL(hydro, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

				return Of_Ks = stats.NASH_SUTCLIFE_MINIMIZE(log1p(hydro.Ks[1:N_iZ]) , log1p(Kₛ_Model[1:N_iZ]))	
			end  # function: OF_KSMODEL
			# --------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : SEARCHRANGE
		#		Required by BlackBoxOptim
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function SEARCHRANGE(optimKsmodel)
				ParamOpt_Min₂ = copy(optim.ParamOpt_Min)
				ParamOpt_Max₂ = copy(optim.ParamOpt_Max)

			return SearchRange = (collect(zip(Float64.(ParamOpt_Min₂), Float64.(ParamOpt_Max₂))))
			end  # function: SEARCHRANGE


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function X_2_τ(ksmodelτ, optimKsmodel, OptLayer, X)

		for iParam = 1:optimKsmodel.NparamOpt
			Paramₐ = X[iParam]
			
			# Getting the current values of every layer of the hydro parameter of interest
				vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[iParam]))

			# Updating the value of the parameters for the soil wanting to optimize by keeping the values constant
				vectParam[OptLayer] = Paramₐ

			# Putting the updated hydro into hydro
				setfield!(ksmodelτ, Symbol(optimKsmodel.ParamOpt[iParam]), vectParam)
		end # for loop

	return ksmodelτ
	end  # function: PARAM
	
end  # module: startKsModel
# ............................................................