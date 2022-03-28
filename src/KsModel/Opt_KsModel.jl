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
		function START_OPT_KSMODEL(GroupBool_Select, hydro, iGroup_Opt, ipGroup, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param)
				
			# Deriving the feasible range of the τ parameters
				SearchRange = SEARCHRANGE(iGroup_Opt, optimKsmodel)

			# Optimisation algorithme, MaxFuncEvals=1000
				Optimization = BlackBoxOptim.bboptimize(X -> OF_KSMODEL(GroupBool_Select, hydro, iGroup_Opt, ipGroup, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[]); SearchRange=SearchRange, NumDimensions=optimKsmodel.NparamOpt[iGroup_Opt], TraceMode=:silent)

			# Deriving the optimal τ parameters from X
				X = BlackBoxOptim.best_candidate(Optimization)

			# Putting X parameters into τ
				ksmodelτ = X_2_τ(iGroup_Opt, ksmodelτ, optimKsmodel, X)

			# Computing optimal KₛModel
				KₛModel = θψ2KsModel.KSMODEL(GroupBool_Select, hydro, ipGroup, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], iGroup_Opt=iGroup_Opt)

	return KₛModel
		end  # function: START_OPT_KSMODEL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_KSMODEL(GroupBool_Select::Vector{Bool}, hydro, iGroup_Opt, ipGroup, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], Unit="MmH", No_Log_Square⍰="Log", WilMot_Ccc⍰="Ccc", KsMinMax=0.005555556)

			# Deriving the optimal τ parameters from X
				ksmodelτ = X_2_τ(iGroup_Opt, ksmodelτ, optimKsmodel, X)

			#	Compuring Ks model
				KₛModel = θψ2KsModel.KSMODEL(GroupBool_Select, hydro, ipGroup, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

			# Computing the Objective function ========================================================
				Ks_ObsTransformed = fill(0.0::Float64, NiZ)
				Ks_SimTransformed = fill(0.0::Float64, NiZ)
			
			# Determening when Ks < KsMinMax
				if option.ksModel.Of_Split_KsSlowKsFast
					KsSmall_True = fill(false, NiZ)
					KsLarge_True = fill(false, NiZ)

					for iZ=1:NiZ
						if GroupBool_Select[iZ] # If we have selected the data
							if hydro.Ks[iZ] ≥ KsMinMax
								KsSmall_True[iZ] = false
								KsLarge_True[iZ] = true
							else
								KsSmall_True[iZ] = true
								KsLarge_True[iZ] = false
							end # if hydro.Ks[iZ] ≥ KsMinMax
						else
							KsSmall_True[iZ] = false
							KsLarge_True[iZ] = false
						end # if GroupBool_Select[iZ]
					end # for iZ=1:NiZ
				end # if option.ksModel.Of_Split_KsSlowKsFast
				#____________________________________

			Ks_ObsTransformed = fill(0.0::Float64, NiZ)
			Ks_SimTransformed = fill(0.0::Float64, NiZ)
			for iZ=1:NiZ
				if Unit == "MmS"
					Ks_ObsTransformed[iZ] = hydro.Ks[iZ]
					Ks_SimTransformed[iZ] = KₛModel[iZ]

				elseif Unit == "MmH"
					Ks_ObsTransformed[iZ] = cst.MmS_2_MmH .* hydro.Ks[iZ]
					Ks_SimTransformed[iZ] = cst.MmS_2_MmH .* KₛModel[iZ]

				else
					error("Unit = $Unit not available" )
				end
			end # for iZ=1:NiZ

			if option.ksModel.Of_Split_KsSlowKsFast
				Ks_ObsTransformed_Small = Ks_ObsTransformed[KsSmall_True[1:NiZ]]
				Ks_SimTransformed_Small = Ks_SimTransformed[KsSmall_True[1:NiZ]]

				Ks_ObsTransformed_Large = Ks_ObsTransformed[KsLarge_True[1:NiZ]]
				Ks_SimTransformed_Large = Ks_SimTransformed[KsLarge_True[1:NiZ]]
			end

			if No_Log_Square⍰ == "Log"
				if option.ksModel.Of_Split_KsSlowKsFast
					Ks_ObsTransformed_Large = log1p.(Ks_ObsTransformed_Large)
					Ks_SimTransformed_Large = log1p.(Ks_SimTransformed_Large)

				else
					Ks_ObsTransformed = log1p.(Ks_ObsTransformed)
					Ks_SimTransformed = log1p.(Ks_SimTransformed)
				end

			elseif No_Log_Square⍰ == "Square"
				if option.ksModel.Of_Split_KsSlowKsFast
					Ks_ObsTransformed_Large = (Ks_ObsTransformed_Large) .^ 0.5
					Ks_SimTransformed_Large = (Ks_SimTransformed_Large) .^ 0.5

				else
					Ks_ObsTransformed = (Ks_ObsTransformed) .^ 0.5
					Ks_SimTransformed = (Ks_SimTransformed) .^ 0.5
				end
			end


			if WilMot_Ccc⍰ == "Wilmot"
				if option.ksModel.Of_Split_KsSlowKsFast
					Of_KsSmall = 1.0 - stats.NSE_WILMOT(Ks_ObsTransformed_Small, Ks_SimTransformed_Small)
					Of_KsLarge = 1.0 - stats.NSE_WILMOT(Ks_ObsTransformed_Large , Ks_SimTransformed_Large)

					Of_Ks = param.ksModel.WeightKsSlow * Of_KsSmall + (1 - param.ksModel.WeightKsSlow) * Of_KsLarge

				else
					Of_Ks = 1.0 - stats.NSE_WILMOT(Ks_ObsTransformed[GroupBool_Select] , Ks_SimTransformed[GroupBool_Select])

				end

			elseif  WilMot_Ccc⍰ == "Ccc"
				if option.ksModel.Of_Split_KsSlowKsFast
					Of_KsSmall = 1.0 - stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(Ks_ObsTransformed_Small , Ks_SimTransformed_Small)

					Of_KsLarge = 1.0 - stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(Ks_ObsTransformed_Large , Ks_SimTransformed_Large)

					Of_Ks = param.ksModel.WeightKsSlow * Of_KsSmall + (1 - param.ksModel.WeightKsSlow) * Of_KsLarge
					
				else
					Of_Ks = 1.0 - stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(Ks_ObsTransformed[GroupBool_Select] , Ks_SimTransformed[GroupBool_Select])
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
		function SEARCHRANGE(iGroup_Opt, optimKsmodel)
			ParamOpt_Min₂ = copy(optimKsmodel.ParamOpt_Min[iGroup_Opt, 1:optimKsmodel.NparamOpt[iGroup_Opt]])
			ParamOpt_Max₂ = copy(optimKsmodel.ParamOpt_Max[iGroup_Opt, 1:optimKsmodel.NparamOpt[iGroup_Opt]])

		return SearchRange = (collect(zip(Float64.(ParamOpt_Min₂), Float64.(ParamOpt_Max₂))))
		end  # function: SEARCHRANGE
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function X_2_τ(iGroup_Opt, ksmodelτ, optimKsmodel, X)
			for iParam = 1:optimKsmodel.NparamOpt[iGroup_Opt]
				Paramₐ = X[iParam]
				
				# Getting the current values of every layer of the hydro parameter of interest
					vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[iGroup_Opt, iParam]))

				# Updating the value of the parameters for the layer wanting to optimize by keeping the other values constant
					vectParam[iGroup_Opt] = Paramₐ

				# Putting the updated hydro into ksmodelτ
					setfield!(ksmodelτ, Symbol(optimKsmodel.ParamOpt[iGroup_Opt, iParam]), vectParam)
			end # for loop
		return ksmodelτ
		end  # function: PARAM
	#..................................................................

	
end  # module: optKsModel
# ========================================================================