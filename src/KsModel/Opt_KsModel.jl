# =============================================================
#		module: optKsModel
# =============================================================
module optKsModel
	import ..stats, ..θψ2KsModel, ..cst
	import BlackBoxOptim
	export START_OPT_KSMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_OPT_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_OPT_KSMODEL(hydro, ipLayer, KₛModel, KₛModel⍰, ksmodelτ, N_iZ, optim, optimKsmodel)
				
			# Deriving the feasible range of the τ parameters
			SearchRange = SEARCHRANGE(ipLayer, optimKsmodel)

			# Optimisation algorithme
				Optimization = BlackBoxOptim.bboptimize(X -> OF_KSMODEL(hydro, ipLayer, KₛModel, KₛModel⍰, ksmodelτ, N_iZ, optim, optimKsmodel, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[]); SearchRange=SearchRange, NumDimensions=optimKsmodel.NparamOpt[ipLayer], TraceMode=:silent, MaxFuncEvals=3000)

				# Deriving the optimal τ parameters from X
					X = BlackBoxOptim.best_candidate(Optimization)

				# Putting X parameters into τ
					ksmodelτ = X_2_τ(ipLayer, ksmodelτ, optimKsmodel, X)

				# Computing optimal KₛModel
					KₛModel = θψ2KsModel.KSMODEL(hydro, KₛModel, KₛModel⍰, ksmodelτ, N_iZ, optim, optimKsmodel; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

		return KₛModel
		end  # function: START_OPT_KSMODEL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_KSMODEL(hydro, ipLayer, KₛModel, KₛModel⍰, ksmodelτ, N_iZ, optim, optimKsmodel, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], Unit="MmH", No_Log_Square⍰="Log", WilMot_Ccc⍰="Ccc", DataSplit=true, KsMinMax=0.005555556, Wsmall=0.54)

			# Deriving the optimal τ parameters from X
				ksmodelτ = X_2_τ(ipLayer, ksmodelτ, optimKsmodel, X)

			#	Compuring Ks model
				KₛModel = θψ2KsModel.KSMODEL(hydro, KₛModel, KₛModel⍰, ksmodelτ, N_iZ, optim, optimKsmodel; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

			# Computing the Objective function ========================================================
				Ks_ObsTransformed = fill(0.0,N_iZ)
				Ks_SimTransformed = fill(0.0,N_iZ)
			
			# Determening when Ks < KsMinMax
				if DataSplit
					KsSmall_True = fill(false, N_iZ)
					KsLarge_True = fill(false, N_iZ)
					for iZ=1:N_iZ
						if hydro.Ks[iZ] ≥ KsMinMax
							KsSmall_True[iZ] = false
							KsLarge_True[iZ] = true
						else
							KsSmall_True[iZ] = true
							KsLarge_True[iZ] = false
						end # if hydro.Ks[iZ] ≥ KsMinMax
					end # for iZ=1:N_iZ
				end # if DataSplit
				#____________________________________

			
			Ks_ObsTransformed = fill(0.0,N_iZ)
			Ks_SimTransformed = fill(0.0,N_iZ)
			for iZ=1:N_iZ
				if Unit == "MmS"
					Ks_ObsTransformed[iZ] = hydro.Ks[iZ]
					Ks_SimTransformed[iZ] = KₛModel[iZ]

				elseif Unit == "MmH"
					Ks_ObsTransformed[iZ] = cst.MmS_2_MmH .* hydro.Ks[iZ]
					Ks_SimTransformed[iZ] = cst.MmS_2_MmH .* KₛModel[iZ]

				else
					error("Unit = $Unit not available" )
				end
			end

			if DataSplit
				Ks_ObsTransformed_Small = Ks_ObsTransformed[KsSmall_True[1:N_iZ]]
				Ks_ObsTransformed_Large = Ks_ObsTransformed[KsLarge_True[1:N_iZ]]

				Ks_SimTransformed_Small = Ks_SimTransformed[KsSmall_True[1:N_iZ]]
				Ks_SimTransformed_Large = Ks_SimTransformed[KsLarge_True[1:N_iZ]]
			end

			if No_Log_Square⍰ == "Log"
				if DataSplit
					# Ks_ObsTransformed_Small = log1p.(Ks_ObsTransformed_Small)
					Ks_ObsTransformed_Large = log1p.(Ks_ObsTransformed_Large)
					# Ks_SimTransformed_Small = log1p.(Ks_SimTransformed_Small)
					Ks_SimTransformed_Large = log1p.(Ks_SimTransformed_Large)

				else
					Ks_ObsTransformed = log1p.(Ks_ObsTransformed)
					Ks_SimTransformed = log1p.(Ks_SimTransformed)
				end

			elseif No_Log_Square⍰ == "Square"
				if DataSplit
					# Ks_ObsTransformed_Small = (Ks_ObsTransformed_Small) .^ 0.5
					Ks_ObsTransformed_Large = (Ks_ObsTransformed_Large) .^ 0.5

					# Ks_SimTransformed_Small = (Ks_SimTransformed_Small) .^ 0.5
					Ks_SimTransformed_Large = (Ks_SimTransformed_Large) .^ 0.5

				else
					Ks_ObsTransformed = (Ks_ObsTransformed) .^ 0.5
					Ks_SimTransformed = (Ks_SimTransformed) .^ 0.5
				end
			end

			if WilMot_Ccc⍰ == "Wilmot"
				if DataSplit == true
					Of_KsSmall = 1.0 - stats.NSE_WILMOT(Ks_ObsTransformed_Small , Ks_SimTransformed_Small)
					Of_KsLarge = 1.0 - stats.NSE_WILMOT(Ks_ObsTransformed_Large , Ks_SimTransformed_Large)

					Of_Ks = Wsmall * Of_KsSmall + (1 - Wsmall) * Of_KsLarge
				
				else
					Of_Ks = 1.0 - stats.NSE_WILMOT(Ks_ObsTransformed[1:N_iZ] , Ks_SimTransformed[1:N_iZ])
				end

			elseif  WilMot_Ccc⍰ == "Ccc"
				if DataSplit == true
					Of_KsSmall = 1.0 - stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(Ks_ObsTransformed_Small , Ks_SimTransformed_Small)

					Of_KsLarge = 1.0 - stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(Ks_ObsTransformed_Large , Ks_SimTransformed_Large)

					Of_Ks = Wsmall * Of_KsSmall + (1 - Wsmall) * Of_KsLarge
					
				else
					Of_Ks = 1.0 - stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(Ks_ObsTransformed[1:N_iZ] , Ks_SimTransformed[1:N_iZ])
				end
			else
				error("WilMot_Ccc⍰ == $WilMot_Ccc⍰ not found")
			end	
		return Of_Ks
		end  # function: OF_KSMODELa
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