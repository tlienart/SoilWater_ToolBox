# =============================================================
#		module: optKsModel
# =============================================================
module optKsModel
	import ..stats, ..θψ2KsModel
	import BlackBoxOptim
	export START_OPT_KSMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_OPT_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_OPT_KSMODEL(hydro, ipLayer, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel)
				
			# Deriving the feasible range of the τ parameters
			SearchRange = SEARCHRANGE(ipLayer, optimKsmodel)

			# Optimisation algorithme
				Optimization = BlackBoxOptim.bboptimize(X -> OF_KSMODEL(hydro, ipLayer, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[]); SearchRange=SearchRange, NumDimensions=optimKsmodel.NparamOpt[ipLayer], TraceMode=:silent)
				# MaxFuncEvals=1500

				# Deriving the optimal τ parameters from X
					X = BlackBoxOptim.best_candidate(Optimization)

				# Putting X parameters into τ
					ksmodelτ = X_2_τ(ipLayer, ksmodelτ, optimKsmodel, X)

				# Computing optimal Kₛ_Model
					Kₛ_Model = θψ2KsModel.KSMODEL(hydro, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

		return Kₛ_Model
		end  # function: START_OPT_KSMODEL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_KSMODEL(hydro, ipLayer, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

			# Deriving the optimal τ parameters from X
				ksmodelτ = X_2_τ(ipLayer, ksmodelτ, optimKsmodel, X)

			#	Compuring Ks model
				Kₛ_Model = θψ2KsModel.KSMODEL(hydro, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

			# Computing the Objective function
				# Of_Ks = stats.NASH_SUTCLIFE_MINIMIZE(log1p.(hydro.Ks[1:N_iZ]) , log1p.(Kₛ_Model[1:N_iZ]))

				# Of_Ks = stats.NASH_SUTCLIFE_MINIMIZE(log.(hydro.Ks[1:N_iZ]) , log.(Kₛ_Model[1:N_iZ]))

				Of_Ks = stats.RMSE(log10.(hydro.Ks[1:N_iZ]) , log10.(Kₛ_Model[1:N_iZ]))

				# println("		====  $Of_Ks ====")

		return Of_Ks
		end  # function: OF_KSMODEL
		# --------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SEARCHRANGE
	#		Required by BlackBoxOptim
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SEARCHRANGE(ipLayer, optimKsmodel)
			ParamOpt_Min₂ = copy(optimKsmodel.ParamOpt_Min[ipLayer, 1:optimKsmodel.NparamOpt[ipLayer]])
			ParamOpt_Max₂ = copy(optimKsmodel.ParamOpt_Max[ipLayer, 1:optimKsmodel.NparamOpt[ipLayer]])

		return SearchRange = (collect(zip(Float64.(ParamOpt_Min₂), Float64.(ParamOpt_Max₂))))
		end  # function: SEARCHRANGE
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function X_2_τ(ipLayer, ksmodelτ, optimKsmodel, X)

			for iParam = 1:optimKsmodel.NparamOpt[ipLayer]

				Paramₐ = X[iParam]
				
				# Getting the current values of every layer of the hydro parameter of interest
					vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[ipLayer, iParam]))

				# Updating the value of the parameters for the layer wanting to optimize by keeping the other values constant
					vectParam[ipLayer] = Paramₐ

				# Putting the updated hydro into ksmodelτ
					setfield!(ksmodelτ, Symbol(optimKsmodel.ParamOpt[ipLayer, iParam]), vectParam)

					# println("		", Symbol(optimKsmodel.ParamOpt[ipLayer, iParam]) , "=" ,vectParam)
			end # for loop

		return ksmodelτ
		end  # function: PARAM
	#..................................................................

	
end  # module: optKsModel
# ========================================================================