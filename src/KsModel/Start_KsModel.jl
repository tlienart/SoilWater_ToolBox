# =============================================================
#		module: startKsModel
# =============================================================

include("θψ_2_KsModel.jl")
include("Opt_KsModel.jl")

module startKsModel
	import ..θψ2KsModel, ..optKsModel, ..stats, ..plot
	export START_KSMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_KSMODEL(hydro, option, param, path, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], ipLayer=1, N_Group=2)

			# GROUPING CLASSES
				GroupBool, ipGroup, N_Group = GROUPING_KSMODEL(hydro, N_Group, NiZ, option, param)

			# OPTIMISE KₛModel
			# If there are τ parameters to be optimised
			if sum(optimKsmodel.NparamOpt) ≥ 1

				for iGroup_Opt=1:N_Group
					println("\n       === iGroup_Opt=$iGroup_Opt === \n")

					# Is there parameters to be optimised in this group
					if optimKsmodel.NparamOpt[iGroup_Opt] ≥ 1

						GroupBool_Select = GroupBool[1:NiZ, iGroup_Opt]

						KₛModel = optKsModel.START_OPT_KSMODEL(GroupBool_Select, hydro, iGroup_Opt, ipGroup, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param)

						ksmodelτ = STATISTICS_KSMODEL(hydro, iGroup_Opt, GroupBool_Select, KₛModel, ksmodelτ, optimKsmodel)


						if option.other.Ploting && option.ksModel.Plot_KsModel
							NameSim = "σ_" * string(iGroup_Opt)

							plot.ksmodel.KSMODEL(hydro, KₛModel[GroupBool_Select], hydro.Ks[GroupBool_Select], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[GroupBool_Select], hydro.θsMacMat[GroupBool_Select], hydro.σ[GroupBool_Select])
						end # if option.Plot

					end # if optimKsmodel.NparamOpt[iGroup_Opt] ≥ 1
				end # for iGroup_Opt=1:N_Group

			# RUN KₛModel
			else
				GroupBool_Select = fill(true, NiZ)

				iGroup_Opt = 1

				KₛModel = θψ2KsModel.KSMODEL(GroupBool_Select, hydro, ipGroup, KₛModel, KₛModel⍰, ksmodelτ, NiZ::Int64, optim, optimKsmodel, option, param; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

				~ = STATISTICS_KSMODEL(hydro, iGroup_Opt, GroupBool_Select, KₛModel, ksmodelτ, optimKsmodel)
			end  # if: optimKsmodel

			if option.other.Ploting
				println("\n       === ALL SIMULATIONS === \n")
				NameSim = "All_"
				plot.ksmodel.KSMODEL(hydro, KₛModel[1:NiZ], hydro.Ks[1:NiZ], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[1:NiZ], hydro.θsMacMat[1:NiZ], hydro.σ[1:NiZ])	
			end
			
			for iZ=1:NiZ
				if "Ks" ∉ optim.ParamOpt
					hydro.Ks[iZ] = KₛModel[iZ]

					# If wanting to assure that the feasible range is physical
					hydro.Ks[iZ] = max( min(hydro.Ks[iZ], hydro.Ks_Max[iZ]), hydro.Ks_Min[iZ])
				end #  hydro.Ks[iZ] < eps(100.0)
			end # if: hydro.Ks[iZ] > eps(10.0)
		
		return hydro, KₛModel
		end  # function: START_KSMODEL
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : GROUPING_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function GROUPING_KSMODEL(hydro, N_Group, NiZ, option, param)
			
         ipGroup   = fill(1::Int64, NiZ)
         GroupBool = fill(true, NiZ, N_Group)

			# Selecting groups if required
			if option.ksModel.Group
				for iZ=1:NiZ
					if hydro.σ[iZ] ≤ param.ksModel.σₛₚₗᵢₜ
						ipGroup[iZ] = 1
						GroupBool[iZ, 1] = true
						GroupBool[iZ, 2] = false
					else
						ipGroup[iZ] = 2
						GroupBool[iZ, 1] = false
						GroupBool[iZ, 2] = true
					end  # if: hydro.  
				end # for iZ=1:NiZ
				N_Group = 2
			else
				N_Group = 1
			end #  option.ksModel.Group
				
		return GroupBool, ipGroup, N_Group
		end  # function: SELECTION
	# ------------------------------------------------------------------
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : STATISTICS_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function STATISTICS_KSMODEL(hydro, iGroup_Opt, GroupBool_Select, KₛModel, ksmodelτ, optimKsmodel)
			# STATISTICS
				ksmodelτ.Nse_τ[iGroup_Opt]    = stats.NSE(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))
				ksmodelτ.Rmse_τ[iGroup_Opt]   = stats.RMSE(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))
				ksmodelτ.Wilmot_τ[iGroup_Opt] = stats.NSE_WILMOT(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))
				ksmodelτ.Ccc_τ[iGroup_Opt]    = stats.stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))

			# PRINING RESULTS
			if sum(optimKsmodel.NparamOpt) ≥ 1 
				println("		 Nse_τ    =  $(ksmodelτ.Nse_τ)")
				println("		 Rmse_τ   =  $(ksmodelτ.Rmse_τ)")
				println("		 Wilmot_τ =  $(ksmodelτ.Wilmot_τ)")
				println("		 Ccc_τ    =  $(ksmodelτ.Ccc_τ)")
			end

				for iParam = 1:optimKsmodel.NparamOpt[iGroup_Opt]	
					# Getting the current values of every layer of the hydro parameter of interest
						vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[iGroup_Opt, iParam]))
						println("		", Symbol(optimKsmodel.ParamOpt[iGroup_Opt, iParam]) , "=" ,vectParam)
				end # for loop
			
		return ksmodelτ
		end  # function: STATISTICS_KSMODEL
	# ------------------------------------------------------------------


end  # module: startKsModel
# =====================================================================