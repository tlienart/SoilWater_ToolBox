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

	# OPTIONS / PARAM / path ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		option = options.OPTIONS()
		param = params.PARAM()
		path = paths.PATH()

	# READING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if option.run.ChangeHydroModel
			# Creating 
			hydroTranslate = hydroStruct.HYDROSTRUCT(1000)
			
			hydroTranslate, N_SoilSelect = reading.READ_STRUCT(hydroTranslate, path.inputSoilwater.ConvertModel)
			
			# Temporary Id
				Id_Select = collect(1:1:N_SoilSelect)
		
			# Deriving a table of θ(Ψ)
				table.hydroLab.TABLE_EXTRAPOINTS_θΨ(hydroTranslate, Id_Select, N_SoilSelect, path.inputSoilwater.Ψθ, param.hydro.Ψ_TableComplete; Orientation="Vertical")
			
			# Deriving a table of K(θ)
				table.hydroLab.TABLE_EXTRAPOINTS_Kθ(hydroTranslate, Id_Select, param.hydro.Ψ_TableComplete, hydroTranslate.Ks[1:N_SoilSelect], N_SoilSelect::Int64, path.inputSoilwater.Kunsat)

			# Creating an Id output required by the program
				table.TABLE_ID(N_SoilSelect::Int64, path.inputSoilwater.Id_Select)
			
		elseif !(option.run.Hypix)
			# Selecting soils of interest
				Id_Select, Id_Select_True, N_SoilSelect = reading.ID(PathIdSlect=path.inputSoilwater.Id_Select, PathOptionSelect=path.option.Select)

			# Determine the soils to simulale
				N_SoilSelect = Int(min(N_SoilSelect, param.globalparam.N_iZ_Simulations))
				Id_Select = Id_Select[1:N_SoilSelect]

			# Reading bulk density
				if option.run.ρb_2_Φ
					RockW, ρ_Rock, ρbSoil, ρp_Fine = reading.BULKDENSITY(Id_Select, N_SoilSelect, path.inputSoilwater.BulkDensity)
				end

			# Reading θ(Ψ)
				if option.run.HydroLabθΨ ≠ :No # <> = <> = <> = <> = <> = <>
					θ_θΨ, Ψ_θΨ, N_θΨ = reading.θΨ(Id_Select, N_SoilSelect, path.inputSoilwater.Ψθ)
				else
					θ_θΨ = [] 
					Ψ_θΨ = [] 
					N_θΨ = 0.0
				end
		
			# Reading K(θ)
				if option.hydro.KunsatΨ # <>=<>=<>=<>=<>
					K_KΨ, Ψ_KΨ, N_KΨ = reading.KUNSATΨ(Id_Select, N_SoilSelect, path)
				else
					Ψ_KΨ = []
					K_KΨ = []
					N_KΨ = 0
				end # option.hydro.KunsatΨ

			# Reading Particle Size distribution
				if option.run.IntergranularMixingPsd
					Rpart, ∑Psd, N_Psd = reading.PSD(Id_Select, N_SoilSelect, path.inputSoilwater.Psd)
				else
					∑Psd = zeros(Float64, N_SoilSelect, 1)
					Rpart = zeros(Float64, N_SoilSelect, 1)
					N_Psd = zeros(Float64, N_SoilSelect)	
				end # option.run.IntergranularMixingPsd
		
			# Reading infiltration	
				if option.run.InfiltBest # <>=<>=<>=<>=<>
					Tinfilt, ∑Infilt_Obs, N_Infilt, infiltParam  = reading.INFILTRATION(Id_Select, N_SoilSelect, path.inputSoilwater.Infiltration, path.inputSoilwater.Infiltration_Param)
				end # option.run.InfiltBest

			# Reading Smap data
				if option.dataFrom.Smap
					smap = readSmap.SMAP(Id_Select_True, N_SoilSelect)
					rfWetable = readSmap.ROCKFRAGMENT_WETTABLE()
				end # option.dataFrom.Smap

			else # TODO: Needs to be removed
				N_SoilSelect = 1
				
			end # Option
			

		# END READING ................................................................
	if option.dataFrom.Smap
		if option.smap.CorrectStone
			println("=\n  				~~~~ Stone correction ~~~~~ \n")
			θ_θΨ = stoneSmap.STONECORRECTION(N_SoilSelect, N_θΨ, smap, θ_θΨ, Ψ_θΨ)
		end
		if option.smap.CorrectStoneWetability
			θ_θΨ = stoneSmap.STONECORRECTION_WETTABLE(N_SoilSelect, N_θΨ, rfWetable, smap, θ_θΨ, Ψ_θΨ)
		end

		if option.smap.UsePointKosugiBimodal
			N_θΨ, θ_θΨ, Ψ_θΨ = reading.θψ_FILE(N_SoilSelect, θ_θΨ, Ψ_θΨ, N_θΨ, path.tableSoilwater.Table_ExtraPoints_θΨ)
		end
	end

	if option.dataFrom.Jules
		SoilName_2_SiteName,  SiteName_2_θini = jules.START_JULES()
		smap2hypix.SMAP_2_HYPIX(SoilName_2_SiteName,  SiteName_2_θini)
		
	end  # if: option.START_JULES()

	if option.run.HydroLabθΨ ≠ :No # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	println("\n === START: DERIVING HYDRO PARAMETERS  === \n")
	println("         ===== Model Name= $(path.option.Model_Name) =====")
		
		# INITIALIZES HYDRAULIC PARAMETERS STRUCT INDEPENDENTLY OF THE SELECTED MODEL
			hydro = hydroStruct.HYDROSTRUCT(option.hydro, N_SoilSelect)

			hydroOther = hydroStruct.HYDRO_OTHERS(N_SoilSelect)

			hydro, optim = reading.HYDRO_PARAM(option.hydro, hydro, N_SoilSelect, path.inputSoilwater.HydroParam_ThetaH)

			checking.CHECKING(option, option.hydro, optim)

			if option.dataFrom.Smap
				hydroParam, optim = stoneSmap.STONECORRECTION_HYDRO(hydro, N_SoilSelect, optim, smap)
			end

			# plotOther.PLOT_σ_2_θr()
			# plotOther.SE_Ψ_CONSTRAINED()
			# plotOther.σ_ψM_SCEARIO()

		# If the hydraulic parameters were already derived than get the data from file instead or rerunning the model	
		if option.run.HydroLabθΨ == :File
			println("    ~ HydroLab HydroParam reading from file ~")
			hydro = reading.HYDROPARAM(Id_Select, N_SoilSelect, hydro)
		else
			# Total Porosity= Φ
				if option.run.ρb_2_Φ
					hydro.Φ = Φ.ρB_2_Φ(N_SoilSelect, RockW, ρ_Rock, ρbSoil, ρp_Fine)
				end

			if option.hydro.KunsatΨ
				hydro, hydroOther = hydrolabOpt.HYDROLABOPT_START(N_SoilSelect=N_SoilSelect, ∑Psd=∑Psd, θ_θΨ=θ_θΨ, Ψ_θΨ=Ψ_θΨ, N_θΨ=N_θΨ, K_KΨ=K_KΨ, Ψ_KΨ=Ψ_KΨ, N_KΨ=N_KΨ, hydro=hydro, hydroOther=hydroOther, option=option, optionₘ=option.hydro, optim=optim)

			else
				hydro, hydroOther =  hydrolabOpt.HYDROLABOPT_START(N_SoilSelect=N_SoilSelect, ∑Psd=∑Psd, θ_θΨ=θ_θΨ, Ψ_θΨ=Ψ_θΨ, N_θΨ=N_θΨ, hydro=hydro, hydroOther=hydroOther, option=option, optionₘ=option.hydro, optim=optim)
			end # option.hydro.KunsatΨ

			# SPATIAL CASE FOR BROOKS AND COREY
				if option.hydro.HydroModel==:BrooksCorey || option.hydro.HydroModel==:ClappHornberger
					for iZ=1:N_SoilSelect
						hydro.Ψga[iZ] = wrc.GREEN_AMPT(optionₘ, iZ, hydro)
					end
				end #  option.hydro.HydroModel
		end # option.run.HydroLabθΨ == :File


		# Deriving Kunsat from θ(Ψ)
		"""Pollacco, J.A.P., Webb, T., McNeill, S., Hu, W., Carrick, S., Hewitt, A., Lilburne, L., 2017. Saturated hydraulic conductivity model computed from bimodal water retention curves for a range of New Zealand soils. Hydrol. Earth Syst. Sci. 21, 2725–2737. https://doi.org/10.5194/hess-21-2725-2017"""
			KunsatModel_Lab = fill(0.0::Float64, N_SoilSelect)
			if option.hydro.HydroModel==:Kosugi 
			println("	=== START: Dering Ks from lab θ(Ψ) data ===")
				for iZ=1:N_SoilSelect
					if option.dataFrom.Smap
						if smap.IsTopsoil[iZ] == 1
							TopsoilSubsoil="Topsoil"
						else
							TopsoilSubsoil="Subsoil"
						end
					else
						TopsoilSubsoil="Topsoil"
					end
					KunsatModel_Lab[iZ] = kunsat.θΨ_2_KUNSAT(option.hydro, 0.9999, iZ, hydro, 0.0; TopsoilSubsoil=TopsoilSubsoil)

					if option.hydro.KunsatΨ == false
						hydro.Ks[iZ] = KunsatModel_Lab[iZ]
					end
				end # for
			println("		~~~ END: Dering Ks from lab θ(Ψ) data ~~~ \n")
			end

	println("=== END  : DERIVING HYDRO PARAMETERS  === \n")
	else
		hydro = []
	end


	if option.run.IntergranularMixingPsd  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	println("=== START: PSD MODEL  ===")
		# Structure of hydroPsd
			hydroPsd = hydroStruct.HYDROSTRUCT(N_SoilSelect)
			hydroOther_Psd = hydroStruct.HYDRO_OTHERS(N_SoilSelect)
			hydroPsd, optim_Psd = reading.HYDRO_PARAM(hydroPsd, N_SoilSelect, path.inputSoilwater.HydroParam_ThetaH)

		# Total Porosity= Φ
		if option.run.ρb_2_Φ
			hydroPsd.Φ = Φ.ρB_2_Φ(N_SoilSelect, RockW, ρ_Rock, ρbSoil, ρp_Fine)
		end

		# PSD model
			paramPsd, N_Psd, θ_Rpart, Ψ_Rpart, Psd, hydroPsd = psdStart.START_PSD(∑Psd, hydro, hydroPsd, N_Psd, N_SoilSelect, N_θΨ, Rpart, θ_θΨ, Ψ_θΨ)

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
			hydroInfilt.Φ = Φ.ρB_2_Φ(N_SoilSelect, RockW, ρ_Rock, ρbSoil, ρp_Fine)

		# Running infiltration model
			infiltOutput, hydroInfilt, ∑Infilt_3D, ∑Infilt_1D = infiltStart.START_INFILTRATION(∑Infilt_Obs, ∑Psd, hydro, hydroInfilt, Id_Select, infiltParam, N_Infilt, N_SoilSelect, Tinfilt)

	println("=== END  : INFILTRATION  === \n")
	else
		hydroInfilt = []
	end # option.run.InfiltBest

	if option.run.Hypix
		hypixStart.HYPIX_START(option)
	end # option.run.Hypix


	# TABLES OUTPUT ======================================================================================
		if option.run.HydroLabθΨ ≠ :No && option.run.HydroLabθΨ ≠ :File # <>=<>=<>=<>=<>

			if !(option.dataFrom.Smap)
				table.hydroLab.θΨK(hydro, hydroOther, Id_Select[1:N_SoilSelect], KunsatModel_Lab, N_SoilSelect, path.Table_θΨK)
			else
				tableSmap.θΨK(hydro, hydroOther, Id_Select[1:N_SoilSelect], KunsatModel_Lab, N_SoilSelect, smap)

				if option.smap.AddPointKosugiBimodal && option.hydro.HydroModel == :Kosugi && option.hydro.σ_2_Ψm == :Constrained
					# Extra points in θ(Ψ) to reduce none uniqueness
					table.hydroLab.TABLE_EXTRAPOINTS_θΨ(option.hydro, hydro, Id_Select, N_SoilSelect, path.tableSoilwater.Table_ExtraPoints_θΨ, param.hydro.Ψ_Table)
	
					# Extra points required by TopNet
						table.hydroLab.TABLE_EXTRAPOINTS_θΨ(option.hydro, hydro, Id_Select, N_SoilSelect, path.tableSoilwater.Table_KosugiθΨ, param.hydro.smap.Ψ_Table)

						table.hydroLab.TABLE_EXTRAPOINTS_Kθ(option.hydro, hydro, Id_Select, param.hydro.K_Table, KunsatModel_Lab, N_SoilSelect, path.inputSoilwater.Kunsat_Model)
				end

				if option.smap.CombineData
					tableSmap.SMAP(option.hydro, Id_Select, N_SoilSelect, smap)
				end
			end

		end # option.run.HydroLabθΨ 

		if option.run.IntergranularMixingPsd # <>=<>=<>=<>=<>
			table.psd.PSD(Id_Select[1:N_SoilSelect], N_SoilSelect, paramPsd, path.tableSoilwater.Table_Psd)

			if option.psd.HydroParam  && option.psd.HydroParam
				table.psd.θΨK_PSD(hydroPsd, Id_Select, KunsatModel_Psd, N_SoilSelect, path.tableSoilwater.Table_Psd)
			end
			
			if option.psd.Table_Psd_θΨ_θ
				table.psd.PSD_θΨ_θ(Id_Select, N_SoilSelect, hydroPsd, path.tableSoilwater.Table_Psd_θΨ_θ)
			end
		end # option.run.IntergranularMixingPsd

		if option.run.InfiltBest # <>=<>=<>=<>=<>
			table.infilt.HYDRO_INFILT(hydroInfilt, Id_Select, KunsatModel_Infilt, N_SoilSelect, path.tableSoilwater.Table_HydroInfilt)

			table.infilt.INFILT(Id_Select, N_SoilSelect, infiltOutput, path.tableSoilwater.Table_Infilt)
		end # option.run.InfiltBest

	# PRINT OUTPUT ======================================================================================
	if option.other.Ploting && !option.run.Hypix
	println("		=== START: PLOTTING  ===")
	
		if option.smap.Plot_Kunsat  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			plotSmap.PLOT_KUNSAT(hydro, N_SoilSelect, smap; N_Se= 1000)
		end

		if option.run.HydroLabθΨ ≠ :No && option.hydro.Plot_θΨ # <>=<>=<>=<>=<>

			if option.dataFrom.Smap
				plotSmap.makie.HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, Id_Select, hydro, KunsatModel_Lab; smap=smap)

				# plotSmap.HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, Id_Select, hydro, KunsatModel_Lab; N_Se=1000, smap=[])
			else
				plot.lab.HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, Id_Select, hydro, KunsatModel_Lab, path.plotSoilwater.Plot_θΨK, path.option.Model_Name)
			end	
		end # option.run.HydroLabθΨ
		if option.run.IntergranularMixingPsd && option.psd.Plot_θr # <>=<>=<>=<>=<>
			plot.psd.PLOT_θr(∑Psd, N_SoilSelect, hydro, hydroPsd, path.plotSoilwater.Plot_Psd_θr, path.plotSoilwater.Plot_IMP_model)
		end
		if option.run.IntergranularMixingPsd && option.psd.Plot_IMP_Model # <>=<>=<>=<>=<>
			plot.psd.PLOT_IMP_MODEL(Id_Select, Rpart, N_Psd, ∑Psd, Psd, N_SoilSelect, hydro, paramPsd) 
		end
		if  option.run.IntergranularMixingPsd && option.psd.Plot_Psd_θΨ && !option.psd.HydroParam
			println("			~ PSD WARNING Sorry cannot plot Plot_Psd_θΨ as option.psd.HydroParam==false ~")
		end
		if option.run.IntergranularMixingPsd && option.psd.Plot_Psd_θΨ && option.psd.HydroParam # <>=<>=<>=<>=<>
			plot.psd.PLOT_PSD_θΨ(Ψ_θΨ, Ψ_Rpart, θ_θΨ, θ_Rpart, N_θΨ, N_SoilSelect, N_Psd, Id_Select, hydroPsd, hydro, path.plotSoilwater.Plot_Psd_θΨ)
		end
		if option.run.InfiltBest && option.infilt.Plot_∑Infilt  # <>=<>=<>=<>=<>
			plot.infilt.PLOT_∑INFILT(Id_Select, N_Infilt, N_SoilSelect, ∑Infilt_Obs, Tinfilt, ∑Infilt_3D, ∑Infilt_1D, infiltOutput, path.plotSoilwater.Plot_∑infilt_Opt )
		end
		# if option.run.InfiltBest && option.infilt.Plot_SeIni_Range # <>=<>=<>=<>=<>
		# Removing GRUtils software to avoid conflict
		# 	# plot.infilt.PLOT_∑INFILT_SEINI(hydroInfilt, Id_Select, infiltOutput, infiltParam, N_SoilSelect)
		# end

		if option.run.InfiltBest && option.infilt.Plot_θΨ
			if option.run.HydroLabθΨ ≠ :No
				plot.infilt.PLOT_∑INFILT_θΨ(hydroInfilt, Id_Select, N_SoilSelect, path.plotSoilwater.Plot_∑infilt_θΨ; hydro=hydro)
			else
				plot.infilt.PLOT_∑INFILT_θΨ(hydroInfilt, Id_Select, N_SoilSelect, path.plotSoilwater.Plot_∑infilt_θΨ)
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