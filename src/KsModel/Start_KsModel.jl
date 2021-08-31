# =============================================================
#		module: startKsModel
# =============================================================
module startKsModel
	import ..θψ2KsModel, ..optKsModel, ..stats
	export START_KSMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_KSMODEL(hydro, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel, optionₘ; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], ipLayer=1)

			# OPTIMISE Kₛ_Model
			if sum(optimKsmodel.NparamOpt) ≥ 1

				for ipLayer = 1:2
					if optimKsmodel.NparamOpt[ipLayer] ≥ 1
						Kₛ_Model = optKsModel.START_OPT_KSMODEL(hydro, ipLayer, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel, optionₘ)
					end
				end

			# RUN Kₛ_Model
			else
				Kₛ_Model = θψ2KsModel.KSMODEL(hydro, Kₛ_Model, ksmodelτ, N_iZ, optim, optimKsmodel, optionₘ; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])
	
			end  # if: optimKsmodel
			
			for iZ=1:N_iZ
				if "Ks" ∉ optim.ParamOpt
					hydro.Ks[iZ] = Kₛ_Model[iZ]
				end #  hydro.Ks[iZ] < eps(100.0)
			end # if: hydro.Ks[iZ] > eps(10.0)


			# STATISTICS
            ksmodelτ.Nse_τ[1]    = stats.NSE(log1p.(hydro.Ks[1:N_iZ]) , log1p.(Kₛ_Model[1:N_iZ]))
            ksmodelτ.Rmse_τ[1]   = stats.RMSE(log1p.(hydro.Ks[1:N_iZ]) , log1p.(Kₛ_Model[1:N_iZ]))
            ksmodelτ.Wilmot_τ[1] = stats.WILMOT(log1p.(hydro.Ks[1:N_iZ]) , log1p.(Kₛ_Model[1:N_iZ]))

				println("		 Nse_τ    =  $(ksmodelτ.Nse_τ)")
				println("		 Rmse_τ   =  $(ksmodelτ.Rmse_τ)")
				println("		 Wilmot_τ =  $(ksmodelτ.Wilmot_τ)")

				for iParam = 1:optimKsmodel.NparamOpt[ipLayer]	
					# Getting the current values of every layer of the hydro parameter of interest
						vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[ipLayer, iParam]))
						println("		", Symbol(optimKsmodel.ParamOpt[ipLayer, iParam]) , "=" ,vectParam)
				end # for loop
				
		return hydro, Kₛ_Model
		end  # function: START_KSMODEL
		#..................................................................
		
	end  # module: startKsModel
	# =====================================================================