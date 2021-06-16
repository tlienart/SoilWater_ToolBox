##========================================================================================
##                                                                                      ##
##                                 Soil Water ToolBox                                   ##
##                                                                                      ##
##========================================================================================

include("Including.jl")

# ===============================================================
#		FUNCTION : START_TOOLBOX
# ==============================================================
function START_TOOLBOX()
	 
	# _______________________ START: option/ param/ path _______________________ 
		
		option = options.OPTIONS()

		param = params.PARAM()

		path = paths.PATH(1, option)

	# ------------------------END: option/ param/ path---------------------------


	# _______________________ START: reading _______________________ 
	println("----- START READING -----------------------------------------------")
		# Determine which soils/ profile to run
			if option.run.Hypix
				IdSelect, IdSelect_True, Soilname, N_SoilSelect = reading.ID(PathIdSlect=path.hyPix.IdSelect, PathOptionSelect=path.option.Select, PathModelName="")
			else
				IdSelect, IdSelect_True, Soilname, N_SoilSelect = reading.ID(PathIdSlect=path.inputSoilwater.IdSelect, PathOptionSelect=path.option.Select, PathModelName=path.option.ModelName)
			end # if: option.run.Hypix

		# If we have bulk density and rock fragment data:
			if option.data.BulkDensity
				RockFragment, ρₚ_Rock, ρᵦ_Soil, ρₚ_Fine = reading.BULKDENSITY(IdSelect, N_SoilSelect, path.inputSoilwater.BulkDensity)
			end

		# if we have Infilt data:
			if option.data.Infilt
				Tinfilt, ∑Infilt_Obs, N_Infilt, infiltParam = reading.INFILTRATION(IdSelect, N_SoilSelect, path.inputSoilwater.Infiltration, path.inputSoilwater.Infiltration_Param)
			end  # if: option.data.Infilt
	
		# If we have θ(Ψ) data:
			if option.data.θΨ
					θ_θΨ, Ψ_θΨ, N_θΨ = reading.θΨ(IdSelect, N_SoilSelect, path.inputSoilwater.Ψθ)
			else
				θ_θΨ = [] 
				Ψ_θΨ = [] 
				N_θΨ = 0	
			end  # if: option.data.θΨ

		# If we have K(θ) data:
			if option.data.Kθ
				K_KΨ, Ψ_KΨ, N_KΨ = reading.KUNSATΨ(IdSelect, N_SoilSelect, path)
			else
				Ψ_KΨ = []
				K_KΨ = []
				N_KΨ = 0
			end  # if: Kθ

		# If we have PSD data:
			if option.data.Psd
				Rpart, ∑Psd, N_Psd = reading.PSD(IdSelect, N_SoilSelect, path.inputSoilwater.Psd)
			else
				∑Psd = []
				# ∑Psd = zeros(Float64, N_SoilSelect, 1)
				Rpart = []
				# Rpart = zeros(Float64, N_SoilSelect, 1)
				# N_Psd = zeros(Float64, N_SoilSelect)	
				N_Psd = 0		
			end  # if: option.data.Psd

		# If we have total porosity data:
			if option.data.TotalPorosity
				RockFragment, Φ = reading.TOTAL_POROSITY(IdSelect, N_SoilSelect, path.inputSoilwater.Φ)	
			end  # if: option.data.TotalPorosity

		# If we have SoilInformation:
			if option.data.SoilInformation
				IsTopsoil, RockClass = reading.SOIL_INFORMATION(IdSelect, N_SoilSelect, path.inputSoilwater.SoilInformation)
			end  # if: option.data.SoilInformation

		#--- NON CORE ----
			# SMAP if we have information of the wetability of rocks:
				if option.data.RockWetability
					rfWetable = reading.smap.ROCKFRAGMENT_WETTABLE(path.inputSmap.LookupTable_RockWetability)	
				end  # if: option.data.RockWetability

			# SMAP if we have information of Smap:
				# if option.data.Smap
				# 	smap = read.smap.SMAP(IdSelect_True, N_SoilSelect, path.inputSmap.Smap)
				# end # option.dataFrom.Smap

	println("----- END READING ----------------------------------------------- \n")
	# ------------------------END: reading---------------------------
	
	# _______________________ START: running HydroLabθΨ _______________________ 
	if option.run.HydroLabθΨ ≠ :No
	println("----- START RUNNING HYDROLABΘΨ -----------------------------------------------")

		# Structures
			hydro = hydroStruct.HYDROSTRUCT(option.hydro, N_SoilSelect)
			hydroOther = hydroStruct.HYDRO_OTHERS(N_SoilSelect)

		# Reading the physical feasible parameter space of the hydraulic parameters
			hydro, optim = reading.HYDRO_PARAM(option.hydro, hydro, N_SoilSelect, path.inputSoilwater.HydroParam_ThetaH)

		# Checking the data
			checking.CHECKING(option, option.hydro, optim)

		# Need to correct the θ(Ψ) curve for rock fragment
		if option.run.RockCorection && option.rockFragment.RockInjectedIncluded==:Injected
			θ_θΨ = rockFragment.injectRock.CORECTION_θΨ!(N_SoilSelect, N_θΨ, RockFragment, θ_θΨ)
		end # if: option.rockFragment.RockInjectedIncluded == :Injected

		# Computing Total Porosity Φ
		if option.run.ρᵦ_2_Φ
			hydro.Φ = rockFragment.ρᵦ_2_Φ(N_SoilSelect, option, RockFragment, ρₚ_Rock, ρᵦ_Soil, ρₚ_Fine)

		elseif option.rockFragment.RockInjectedIncluded == :Injected && option.data.TotalPorosity
			hydro.Φ = rockFragment.injectRock.CORECTION_Φ!(N_SoilSelect, RockFragment, Φ)
		end

		if option.hydro.KunsatΨ
			hydro, hydroOther = hydrolabOpt.HYDROLABOPT_START(N_SoilSelect=N_SoilSelect, ∑Psd=∑Psd, θ_θΨ=θ_θΨ, Ψ_θΨ=Ψ_θΨ, N_θΨ=N_θΨ, K_KΨ=K_KΨ, Ψ_KΨ=Ψ_KΨ, N_KΨ=N_KΨ, hydro=hydro, hydroOther=hydroOther, option=option, optionₘ=option.hydro, optim=optim, param=param)

		else
			hydro, hydroOther = hydrolabOpt.HYDROLABOPT_START(N_SoilSelect=N_SoilSelect, ∑Psd=∑Psd, θ_θΨ=θ_θΨ, Ψ_θΨ=Ψ_θΨ, N_θΨ=N_θΨ, hydro=hydro, hydroOther=hydroOther, option=option, optionₘ=option.hydro, optim=optim, param=param)
		end # option.hydro.KunsatΨ

		# Spatial case
		if option.hydro.HydroModel==:BrooksCorey || option.hydro.HydroModel==:ClappHornberger
			for iZ=1:N_SoilSelect
				hydro.Ψga[iZ] = wrc.GREEN_AMPT(option.hydro, iZ, hydro)
			end
		end #  option.hydro.HydroModel

	println("----- END: RUNNING HYDROLABΘΨ ----------------------------------------------- \n")
	else
			hydro = []
	end # option.run.HydroLabθΨ

# If the hydraulic parameters were already derived than get the data from file instead or rerunning the model	
		# if option.run.HydroLabθΨ == :File
		# 	println("    ~ HydroLab HydroParam reading from file ~")
		# 	hydro = reading.HYDROPARAM(IdSelect, N_SoilSelect, hydro)
		# else
				# if option.dataFrom.Smap
		# 	hydroParam, optim = stoneSmap.STONECORRECTION_HYDRO(hydro, N_SoilSelect, optim, smap)
		# end
		# ------------------------END: running HydroLab---------------------------  
	
	

		if option.run.ChangeHydroModel
			# Creating 
			hydroTranslate = hydroStruct.HYDROSTRUCT(1000)
			
			hydroTranslate, N_SoilSelect = reading.READ_STRUCT(hydroTranslate, path.inputSoilwater.ConvertModel)
			
			# Temporary Id
				IdSelect = collect(1:1:N_SoilSelect)
		
			# Deriving a table of θ(Ψ)
				table.hydroLab.TABLE_EXTRAPOINTS_θΨ(hydroTranslate, IdSelect, N_SoilSelect, path.inputSoilwater.Ψθ, param.hydro.Ψ_TableComplete; Orientation="Vertical")
			
			# Deriving a table of K(θ)
				table.hydroLab.TABLE_EXTRAPOINTS_Kθ(hydroTranslate, IdSelect, param.hydro.Ψ_TableComplete, hydroTranslate.Ks[1:N_SoilSelect], N_SoilSelect::Int64, path.inputSoilwater.Kunsat)

			# Creating an Id output required by the program
				table.TABLE_ID(N_SoilSelect::Int64, path.inputSoilwater.IdSelect)
			
		elseif !(option.run.Hypix)
			else # TODO: Needs to be removed
				N_SoilSelect = 1
			end # Option
			


	if option.dataFrom.Smap
		if option.smap.CorrectStone
			println("=\n  				~~~~ Stone correction ~~~~~ \n")
			θ_θΨ = rockFragment.injectRock.CORECTION_θΨ!(N_SoilSelect, N_θΨ, smap.RockFragment, θ_θΨ)
		end
		if option.smap.CorrectStoneWetability
			θ_θΨ = stoneSmap.STONECORRECTION_WETTABLE(N_SoilSelect, N_θΨ, rfWetable, smap, θ_θΨ, Ψ_θΨ)
		end

		if option.smap.UsePointKosugiBimodal
			N_θΨ, θ_θΨ, Ψ_θΨ = reading.θψ_FILE(
				N_SoilSelect, N_θΨ, param, path.tableSoilwater.Table_ExtraPoints_θΨ, θ_θΨ, Ψ_θΨ)
		end
	end

	if option.dataFrom.Jules
		SoilName_2_SiteName,  SiteName_2_θini = jules.START_JULES(path)
		smap2hypix.SMAP_2_HYPIX(SoilName_2_SiteName, SiteName_2_θini, path)	
	end  # if: option.START_JULES()

	

		# Deriving Kunsat from θ(Ψ)
		# """Pollacco, J.A.P., Webb, T., McNeill, S., Hu, W., Carrick, S., Hewitt, A., Lilburne, L., 2017. Saturated hydraulic conductivity model computed from bimodal water retention curves for a range of New Zealand soils. Hydrol. Earth Syst. Sci. 21, 2725–2737. https://doi.org/10.5194/hess-21-2725-2017"""

		# 	if option.hydro.HydroModel == :Kosugi
		# 	RockFragment = 0.0
		# 		for iZ=1:N_SoilSelect
		# 			if hydro.Ks[iZ] > eps(10.0)
		# 				if option.dataFrom.Smap
		# 					hydro.Ks[iZ] = θψ2Ks.θΨ_2_KS(hydro, iZ, option, param, RockFragment[iZ]; IsTopsoil=smap.IsTopsoil[iZ])
		# 				else
		# 					hydro.Ks[iZ] = θψ2Ks.θΨ_2_KS(hydro, iZ, option, param, RockFragment[iZ]; IsTopsoil=1)
		# 				end # if: option.dataFrom.Smap
		# 			end # if: hydro.Ks[iZ] > eps(10.0)
		# 		end # for; iZ=1:N_SoilSelect
		# 	end # if:  option.hydro.HydroModel == :Kosugi 





	if option.run.IntergranularMixingPsd  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	println("=== START: PSD MODEL  ===")
		# Structure of hydroPsd
			hydroPsd = hydroStruct.HYDROSTRUCT(N_SoilSelect)
			hydroOther_Psd = hydroStruct.HYDRO_OTHERS(N_SoilSelect)
			hydroPsd, optim_Psd = reading.HYDRO_PARAM(hydroPsd, N_SoilSelect, path.inputSoilwater.HydroParam_ThetaH)

		# Total Porosity= Φ
		if option.run.ρᵦ_2_Φ
			hydroPsd.Φ = rockFragment.ρᵦ_2_Φ(N_SoilSelect, option, RockFragment, ρₚ_Rock, ρᵦ_Soil, ρₚ_Fine)
		end

		# PSD model
			paramPsd, N_Psd, θ_Rpart, Ψ_Rpart, Psd, hydroPsd = psdStart.START_PSD(∑Psd, hydro, hydroPsd, N_Psd, N_SoilSelect, N_θΨ, param, Rpart, θ_θΨ, Ψ_θΨ)

		KunsatModel_Psd = fill(0.0::Float64, N_SoilSelect)

		if  option.psd.HydroParam
			hydroPsd, hydroOther_Psd = hydrolabOpt.HYDROLABOPT_START(N_SoilSelect=N_SoilSelect, ∑Psd=∑Psd, θ_θΨ=θ_Rpart, Ψ_θΨ=Ψ_Rpart, N_θΨ=N_Psd, hydro=hydroPsd, hydroOther=hydroOther_Psd, option=option, optionₘ=option.psd, optim=optim_Psd)
		end

	
		println("=== END  : DERIVING HYDRO PARAMETERS  === \n")

	println("=== END : PSD MODEL  === \n")
	else
		θ_Rpart = zeros(Float64, N_SoilSelect,1)
		Ψ_Rpart = zeros(Float64, N_SoilSelect,1)
		hydroPsd = hydroStruct.HYDROSTRUCT(option.psd, N_SoilSelect)
		N_Psd = zeros(Float64, N_SoilSelect)

	end # option.run.IntergranularMixingPsd ...............................................................................

	
	if option.run.InfiltBest  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	println("=== START: INFILTRATION  ===")
		# Structure of hydroInfilt
			hydroInfilt = hydroStruct.HYDROSTRUCT(N_SoilSelect)
			hydroOther_Infilt = hydroStruct.HYDRO_OTHERS(N_SoilSelect)
			hydroInfilt, optim_Infilt = reading.HYDRO_PARAM(hydroPsd, N_SoilSelect, path.inputSoilwater.HydroParam_Infilt)

		# Total Porosity= Φ
			hydroInfilt.Φ = rockFragment.ρᵦ_2_Φ(N_SoilSelect, option, RockFragment, ρₚ_Rock, ρᵦ_Soil, ρₚ_Fine)

		# Running infiltration model
			infiltOutput, hydroInfilt, ∑Infilt_3D, ∑Infilt_1D = infiltStart.START_INFILTRATION(∑Infilt_Obs, ∑Psd, hydro, hydroInfilt, IdSelect, infiltParam, N_Infilt, N_SoilSelect, Tinfilt)

	println("=== END  : INFILTRATION  === \n")
	else
		hydroInfilt = []
	end # option.run.InfiltBest

	if option.run.Hypix
		hypixStart.HYPIX_START(Soilname, option, param, path)
	end # option.run.Hypix


	# TABLES OUTPUT ======================================================================================
		if option.run.HydroLabθΨ ≠ :No && option.run.HydroLabθΨ ≠ :File # <>=<>=<>=<>=<>
			KunsatModel_Lab = [] # Temporary
			if !(option.dataFrom.Smap)
				table.hydroLab.θΨK(hydro, hydroOther, IdSelect[1:N_SoilSelect], KunsatModel_Lab, N_SoilSelect, path.tableSoilwater.Table_θΨK)
			else
				tableSmap.θΨK(hydro, hydroOther, IdSelect[1:N_SoilSelect], KunsatModel_Lab, N_SoilSelect, smap, path.Table_θΨK)

				if option.smap.AddPointKosugiBimodal && option.hydro.HydroModel == :Kosugi && option.hydro.σ_2_Ψm == :Constrained
					# Extra points in θ(Ψ) to reduce none uniqueness
					table.hydroLab.TABLE_EXTRAPOINTS_θΨ(option.hydro, hydro, IdSelect, N_SoilSelect, path.tableSoilwater.Table_ExtraPoints_θΨ, param.hydro.Ψ_Table)
	
					# Extra points required by TopNet
						table.hydroLab.TABLE_EXTRAPOINTS_θΨ(option.hydro, hydro, IdSelect, N_SoilSelect, path.tableSoilwater.Table_KosugiθΨ, param.hydro.smap.Ψ_Table)

						table.hydroLab.TABLE_EXTRAPOINTS_Kθ(option.hydro, hydro, IdSelect, param.hydro.K_Table, KunsatModel_Lab, N_SoilSelect, path.inputSoilwater.Kunsat_Model)
				end

				if option.smap.CombineData
					tableSmap.SMAP(option.hydro, IdSelect, N_SoilSelect, smap, param, path)
				end
			end

		end # option.run.HydroLabθΨ 

		if option.run.IntergranularMixingPsd # <>=<>=<>=<>=<>
			table.psd.PSD(IdSelect[1:N_SoilSelect], N_SoilSelect, paramPsd, path.tableSoilwater.Table_Psd)

			if option.psd.HydroParam  && option.psd.HydroParam
				table.psd.θΨK_PSD(hydroPsd, IdSelect, KunsatModel_Psd, N_SoilSelect, path.tableSoilwater.Table_Psd)
			end
			
			if option.psd.Table_Psd_θΨ_θ
				table.psd.PSD_θΨ_θ(IdSelect, N_SoilSelect, hydroPsd, param, path.tableSoilwater.Table_Psd_θΨ_θ)
			end
		end # option.run.IntergranularMixingPsd

		if option.run.InfiltBest # <>=<>=<>=<>=<>
			table.infilt.HYDRO_INFILT(hydroInfilt, IdSelect, KunsatModel_Infilt, N_SoilSelect, path.tableSoilwater.Table_HydroInfilt)

			table.infilt.INFILT(IdSelect, N_SoilSelect, infiltOutput, path.tableSoilwater.Table_Infilt)
		end # option.run.InfiltBest

	# PRINT OUTPUT ======================================================================================
	if option.other.Ploting && !option.run.Hypix
	println("		=== START: PLOTTING  ===")
	
		# if option.smap.Plot_Kunsat  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	plotSmap.PLOT_KUNSAT(hydro, N_SoilSelect, smap; N_Se= 1000)
		# end

		if option.run.HydroLabθΨ ≠ :No && option.hydro.Plot_θΨ # <>=<>=<>=<>=<>

			if option.dataFrom.Smap
				plotSmap.makie.HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, IdSelect, hydro, KunsatModel_Lab, path; smap=smap)

				# plotSmap.HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, IdSelect, hydro, KunsatModel_Lab; N_Se=1000, smap=[])
			else
				plot.lab.HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, IdSelect, hydro, KunsatModel_Lab, path.plotSoilwater.Plot_θΨK, path.option.ModelName)
			end	
		end # option.run.HydroLabθΨ
		if option.run.IntergranularMixingPsd && option.psd.Plot_θr # <>=<>=<>=<>=<>
			plot.psd.PLOT_θr(∑Psd, N_SoilSelect, hydro, hydroPsd, path.plotSoilwater.Plot_Psd_θr, path.plotSoilwater.Plot_IMP_model)
		end
		if option.run.IntergranularMixingPsd && option.psd.Plot_IMP_Model # <>=<>=<>=<>=<>
			plot.psd.PLOT_IMP_MODEL(IdSelect, Rpart, N_Psd, ∑Psd, Psd, N_SoilSelect, hydro, paramPsd) 
		end
		if  option.run.IntergranularMixingPsd && option.psd.Plot_Psd_θΨ && !option.psd.HydroParam
			println("			~ PSD WARNING Sorry cannot plot Plot_Psd_θΨ as option.psd.HydroParam==false ~")
		end
		if option.run.IntergranularMixingPsd && option.psd.Plot_Psd_θΨ && option.psd.HydroParam # <>=<>=<>=<>=<>
			plot.psd.PLOT_PSD_θΨ(Ψ_θΨ, Ψ_Rpart, θ_θΨ, θ_Rpart, N_θΨ, N_SoilSelect, N_Psd, IdSelect, hydroPsd, hydro, path.plotSoilwater.Plot_Psd_θΨ)
		end
		if option.run.InfiltBest && option.infilt.Plot_∑Infilt  # <>=<>=<>=<>=<>
			plot.infilt.PLOT_∑INFILT(IdSelect, N_Infilt, N_SoilSelect, ∑Infilt_Obs, Tinfilt, ∑Infilt_3D, ∑Infilt_1D, infiltOutput, path.plotSoilwater.Plot_∑infilt_Opt )
		end
		# if option.run.InfiltBest && option.infilt.Plot_SeIni_Range # <>=<>=<>=<>=<>
		# Removing GRUtils software to avoid conflict
		# 	# plot.infilt.PLOT_∑INFILT_SEINI(hydroInfilt, IdSelect, infiltOutput, infiltParam, N_SoilSelect)
		# end

		if option.run.InfiltBest && option.infilt.Plot_θΨ
			if option.run.HydroLabθΨ ≠ :No
				plot.infilt.PLOT_∑INFILT_θΨ(hydroInfilt, IdSelect, N_SoilSelect, path.plotSoilwater.Plot_∑infilt_θΨ; hydro=hydro)
			else
				plot.infilt.PLOT_∑INFILT_θΨ(hydroInfilt, IdSelect, N_SoilSelect, path.plotSoilwater.Plot_∑infilt_θΨ)
			end # option.run.HydroLabθΨ
		end # option.run.InfiltBest

	println("=== END: PLOTTING  === \n")
	end # if option.other.Ploting

	# Playing sounds...
		println("\007")

end  # function: START_TOOLBOX
# ..............................................................

println("\n\n ===== START SOIL WATER TOOLBOX =====")
	@time START_TOOLBOX()
println("==== END SOIL WATER TOOLBOX ====")