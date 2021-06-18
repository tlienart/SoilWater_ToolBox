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
				IdSelect, IdSelect_True, Soilname, N_iZ = reading.ID(PathIdSlect=path.hyPix.IdSelect, PathOptionSelect=path.option.Select, PathModelName="")
			else
				IdSelect, IdSelect_True, Soilname, N_iZ = reading.ID(PathIdSlect=path.inputSoilwater.IdSelect, PathOptionSelect=path.option.Select, PathModelName=path.option.ModelName)
			end # if: option.run.Hypix

		# If we have bulk density and rock fragment data:
			if option.data.Φmethod == :ρᵦ
				RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil = reading.BULKDENSITY(IdSelect, N_iZ, path.inputSoilwater.BulkDensity)

				# Φ  corrected for RockFragments
				Φ = rockFragment.ρᵦ_2_Φ(N_iZ, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)
			elseif option.data.Φmethod == :Φ # Total Porosity
				RockFragment, Φ = reading.Φ(IdSelect, N_iZ, path.inputSoilwater.Φ)
				
				Φ = rockFragment.injectRock.CORECTION_Φ!(N_iZ, option, RockFragment, Φ)	
			end # option.data.Φmethod == :ρᵦ

		# if we have Infilt data:
			if option.data.Infilt
				Tinfilt, ∑Infilt_Obs, N_Infilt, infiltParam = reading.INFILTRATION(IdSelect, N_iZ, path.inputSoilwater.Infiltration, path.inputSoilwater.Infiltration_Param)
			end  # if: option.data.Infilt
	
		# If we have θ(Ψ) data:
			if option.data.θΨ
					θ_θΨobs, Ψ_θΨobs, N_θΨobs = reading.θΨ(IdSelect, N_iZ, path.inputSoilwater.Ψθ)
			end  # if: option.data.θΨ

		# If we have K(θ) data:
			if option.data.Kθ
				K_KΨobs, Ψ_KΨobs, N_KΨobs = reading.KUNSATΨ(IdSelect, N_iZ, path)
			end  # if: Kθ

		# If we have PSD data:
			if option.data.Psd
				Rpart, ∑Psd, N_Psd = reading.PSD(IdSelect, N_iZ, path.inputSoilwater.Psd)
			else
				∑Psd = []
				# ∑Psd = zeros(Float64, N_iZ, 1)
				Rpart = []
				# Rpart = zeros(Float64, N_iZ, 1)
				# N_Psd = zeros(Float64, N_iZ)	
				N_Psd = 0		
			end  # if: option.data.Psd

		# If we have SoilInformation:
			if option.data.SoilInformation
				IsTopsoil, RockClass = reading.SOIL_INFORMATION(IdSelect, N_iZ, path.inputSoilwater.SoilInformation)
			end  # if: option.data.SoilInformation

		#--- NON CORE ----
			# SMAP if we have information of the wetability of rocks:
				if option.data.RockWetability
					rfWetable = reading.smap.ROCKFRAGMENT_WETTABLE(path.inputSmap.LookupTable_RockWetability)	
				end  # if: option.data.RockWetability

			# SMAP if we have information of Smap:
				# if option.data.Smap
				# 	smap = read.smap.SMAP(IdSelect_True, N_iZ, path.inputSmap.Smap)
				# end # option.dataFrom.Smap

	println("----- END READING ----------------------------------------------- \n")
	# ------------------------END: reading---------------------------
	
	# _______________________ START: running HydroLabθΨ _______________________ 
	if option.run.HydroLabθΨ ≠ :No
	println("----- START RUNNING HYDROLABΘΨ -----------------------------------------------")

		# STRUCTURES
			hydro = hydroStruct.HYDROSTRUCT(option.hydro, N_iZ)
			hydroOther = hydroStruct.HYDRO_OTHERS(N_iZ)
			hydro, optim = reading.HYDRO_PARAM(option.hydro, hydro, N_iZ, path.inputSoilwater.HydroParam_ThetaH)

		# CHECKING THE DATA
			checking.CHECKING(option, option.hydro, optim)

		# TRANSFERING Φ -> hydro
			if option.data.Φmethod ≠ :No
				for iZ =1:N_iZ 
					hydro.Φ[iZ] = Φ[iZ]
				end
			end # option.data.Φmethod ≠ :No

		# CORRECT θ(Ψ) FOR ROCK FRAGMENT
		if option.run.RockCorection && option.rockFragment.RockInjectedIncluded==:InjectRock
			θ_θΨobs = rockFragment.injectRock.CORECTION_θΨ!(N_iZ, N_θΨobs, RockFragment, θ_θΨobs)
		end # if: option.rockFragment.RockInjectedIncluded == :InjectRock

		# OPTIMISING THE HYDRAULIC PARAMETERS
		if option.hydro.KsOpt
			hydro, hydroOther = hydrolabOpt.HYDROLABOPT_START(N_iZ=N_iZ, ∑Psd=∑Psd, θ_θΨobs=θ_θΨobs, Ψ_θΨobs=Ψ_θΨobs, N_θΨobs=N_θΨobs, K_KΨobs=K_KΨobs, Ψ_KΨobs=Ψ_KΨobs, N_KΨobs=N_KΨobs, hydro=hydro, hydroOther=hydroOther, option=option, optionₘ=option.hydro, optim=optim, param=param)

		else
			hydro, hydroOther = hydrolabOpt.HYDROLABOPT_START(N_iZ=N_iZ, ∑Psd=∑Psd, θ_θΨobs=θ_θΨobs, Ψ_θΨobs=Ψ_θΨobs, N_θΨobs=N_θΨobs, hydro=hydro, hydroOther=hydroOther, option=option, optionₘ=option.hydro, optim=optim, param=param)
		end # option.hydro.KsOpt

		# COMPUTE KS FROM Θ(Ψ)
			Kₛ_Model = fill(0.0::Float64, N_iZ)
			if option.hydro.HydroModel == :Kosugi
				println("\n	=== === Computing model Ks === === ")

				IsTopsoil₁ = 0.0 ; RockFragment₁ = 0.0
				for iZ=1:N_iZ
					if @isdefined RockFragment
						RockFragment₁ = RockFragment[iZ]
					end #@isdefined RockFragment
					if @isdefined IsTopsoil
						IsTopsoil₁ = IsTopsoil[iZ]					
					end  # if: @isdefined IsTopsoil

					Kₛ_Model[iZ] = θψ2Ks.θΨ_2_KS(hydro, iZ, param; RockFragment=RockFragment₁, IsTopsoil=IsTopsoil₁)

					if hydro.Ks[iZ] < eps(100.0)
						hydro.Ks[iZ] = Kₛ_Model[iZ]
					end #  hydro.Ks[iZ] < eps(100.0)
				end # if: hydro.Ks[iZ] > eps(10.0)
				println("	=== === ~~~~~~~~~~~~~~~~~~~~~~~~ === === ")
			end # if: option.hydro.HydroModel == :Kosugi 


		# SPECIAL CASE
			if option.hydro.HydroModel==:BrooksCorey || option.hydro.HydroModel==:ClappHornberger
				for iZ=1:N_iZ
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
		# 	hydro = reading.HYDROPARAM(IdSelect, N_iZ, hydro)
		# else
				# if option.dataFrom.Smap

		# end
		# ------------------------END: running HydroLab---------------------------  
	
	

		if option.run.ChangeHydroModel
			# Creating 
			hydroTranslate = hydroStruct.HYDROSTRUCT(1000)
			
			hydroTranslate, N_iZ = reading.READ_STRUCT(hydroTranslate, path.inputSoilwater.ConvertModel)
			
			# Temporary Id
				IdSelect = collect(1:1:N_iZ)
		
			# Deriving a table of θ(Ψ)
				table.hydroLab.TABLE_EXTRAPOINTS_θΨ(hydroTranslate, IdSelect, N_iZ, path.inputSoilwater.Ψθ, param.hydro.TableComplete_θΨ; Orientation="Vertical")
			
			# Deriving a table of K(θ)
				table.hydroLab.TABLE_EXTRAPOINTS_Kθ(hydroTranslate, IdSelect, param.hydro.TableComplete_θΨ, hydroTranslate.Ks[1:N_iZ], N_iZ::Int64, path.inputSoilwater.Kunsat)

			# Creating an Id output required by the program
				table.TABLE_ID(N_iZ::Int64, path.inputSoilwater.IdSelect)
			
		elseif !(option.run.Hypix)
			else # TODO: Needs to be removed
				N_iZ = 1
			end # Option
			

	if option.dataFrom.Smap
		if option.smap.CorrectStone
			println("=\n  				~~~~ Stone correction ~~~~~ \n")
			θ_θΨobs = rockFragment.injectRock.CORECTION_θΨ!(N_iZ, N_θΨobs, smap.RockFragment, θ_θΨobs)
		end
		if option.smap.CorrectStoneWetability
			θ_θΨobs = stoneSmap.CORECTION_θΨ_WETABLE!(N_iZ, N_θΨobs, rfWetable, smap, θ_θΨobs, Ψ_θΨobs)
		end

		if option.smap.UsePointKosugiBimodal
			N_θΨobs, θ_θΨobs, Ψ_θΨobs = reading.θψ_FILE(
				N_iZ, N_θΨobs, param, path.tableSoilwater.TableExtraPoints_θΨ, θ_θΨobs, Ψ_θΨobs)
		end
	end

	if option.dataFrom.Jules
		SoilName_2_SiteName,  SiteName_2_θini = jules.START_JULES(path)
		smap2hypix.SMAP_2_HYPIX(SoilName_2_SiteName, SiteName_2_θini, path)	
	end  # if: option.START_JULES()

	


	if option.run.IntergranularMixingPsd  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	println("=== START: PSD MODEL  ===")
		# Structure of hydroPsd
			hydroPsd = hydroStruct.HYDROSTRUCT(N_iZ)
			hydroOther_Psd = hydroStruct.HYDRO_OTHERS(N_iZ)
			hydroPsd, optim_Psd = reading.HYDRO_PARAM(hydroPsd, N_iZ, path.inputSoilwater.HydroParam_ThetaH)

		# Total Porosity= Φ
		# if option.run.ρᵦ_2_Φ
			hydroPsd.Φ = rockFragment.ρᵦ_2_Φ(N_iZ, option, RockFragment, ρₚ_Rock, ρᵦ_Soil, ρₚ_Fine)
		# end

		# PSD model
			paramPsd, N_Psd, θ_Rpart, Ψ_Rpart, Psd, hydroPsd = psdStart.START_PSD(∑Psd, hydro, hydroPsd, N_Psd, N_iZ, N_θΨobs, param, Rpart, θ_θΨobs, Ψ_θΨobs)

		KunsatModel_Psd = fill(0.0::Float64, N_iZ)

		if  option.psd.HydroParam
			hydroPsd, hydroOther_Psd = hydrolabOpt.HYDROLABOPT_START(N_iZ=N_iZ, ∑Psd=∑Psd, θ_θΨobs=θ_Rpart, Ψ_θΨobs=Ψ_Rpart, N_θΨobs=N_Psd, hydro=hydroPsd, hydroOther=hydroOther_Psd, option=option, optionₘ=option.psd, optim=optim_Psd)
		end

	
		println("=== END  : DERIVING HYDRO PARAMETERS  === \n")

	println("=== END : PSD MODEL  === \n")
	else
		θ_Rpart = zeros(Float64, N_iZ,1)
		Ψ_Rpart = zeros(Float64, N_iZ,1)
		hydroPsd = hydroStruct.HYDROSTRUCT(option.psd, N_iZ)
		N_Psd = zeros(Float64, N_iZ)

	end # option.run.IntergranularMixingPsd ...............................................................................

	
	if option.run.InfiltBest  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	println("=== START: INFILTRATION  ===")
		# Structure of hydroInfilt
			hydroInfilt = hydroStruct.HYDROSTRUCT(N_iZ)
			hydroOther_Infilt = hydroStruct.HYDRO_OTHERS(N_iZ)
			hydroInfilt, optim_Infilt = reading.HYDRO_PARAM(hydroPsd, N_iZ, path.inputSoilwater.HydroParam_Infilt)

		# Total Porosity= Φ
			hydroInfilt.Φ = rockFragment.ρᵦ_2_Φ(N_iZ, option, RockFragment, ρₚ_Rock, ρᵦ_Soil, ρₚ_Fine)

		# Running infiltration model
			infiltOutput, hydroInfilt, ∑Infilt_3D, ∑Infilt_1D = infiltStart.START_INFILTRATION(∑Infilt_Obs, ∑Psd, hydro, hydroInfilt, IdSelect, infiltParam, N_Infilt, N_iZ, Tinfilt)

	println("=== END  : INFILTRATION  === \n")
	else
		hydroInfilt = []
	end # option.run.InfiltBest

	if option.run.Hypix
		hypixStart.HYPIX_START(Soilname, option, param, path)
	end # option.run.Hypix

	 # _______________________ START: table _______________________ 

	 if option.run.HydroLabθΨ ≠ :No && option.run.HydroLabθΨ ≠ :File # <>=<>=<>=<>=<>
		if !(option.dataFrom.Smap)

			# CORE OUTPUT
				table.hydroLab.θΨK(hydro, hydroOther, IdSelect[1:N_iZ], Kₛ_Model[1:N_iZ], N_iZ, path.tableSoilwater.Table_θΨK)

			# Adding extra points to θ(Ψ) if not available
				table.hydroLab.TABLE_EXTRAPOINTS_θΨ(option.hydro, hydro, IdSelect, N_iZ, path.tableSoilwater.TableExtraPoints_θΨ, param.hydro.TableExtraPoints_θΨ)
			
			# Extra points required by other models
				table.hydroLab.TABLE_EXTRAPOINTS_θΨ(option.hydro, hydro, IdSelect, N_iZ, path.tableSoilwater.TableComplete_θΨ, param.hydro.TableComplete_θΨ; Orientation="Horizontal")

				# table.hydroLab.TABLE_EXTRAPOINTS_Kθ(option.hydro, hydro, IdSelect, param.hydro.K_Table, Kₛ_Model, N_iZ, path.inputSoilwater.Kunsat_Model)

		else
			tableSmap.θΨK(hydro, hydroOther, IdSelect[1:N_iZ], Kₛ_Model, N_iZ, smap, path.Table_θΨK)

			if option.smap.AddPointKosugiBimodal && option.hydro.HydroModel == :Kosugi && option.hydro.σ_2_Ψm == :Constrained
				# Extra points in θ(Ψ) to reduce none uniqueness



			end
			if option.smap.CombineData
				tableSmap.SMAP(option.hydro, IdSelect, N_iZ, smap, param, path)
			end
		end
	 # ------------------------END: table--------------------------- 
	 
	 # _______________________ START: plotting _______________________ 

	 if option.other.Ploting && !option.run.Hypix
	println("		=== START: PLOTTING  ===")
		plot.lab.HYDROPARAM(hydro, IdSelect, K_KΨobs, N_iZ, N_KΨobs, N_θΨobs, option, option.hydro, param, path, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs; N_Se=1000)

	println("		=== END: PLOTTING  === \n")
	 end
	

	 # ------------------------END: plotting---------------------------  




		end # option.run.HydroLabθΨ 

		if option.run.IntergranularMixingPsd # <>=<>=<>=<>=<>
			table.psd.PSD(IdSelect[1:N_iZ], N_iZ, paramPsd, path.tableSoilwater.Table_Psd)

			if option.psd.HydroParam  && option.psd.HydroParam
				table.psd.θΨK_PSD(hydroPsd, IdSelect, KunsatModel_Psd, N_iZ, path.tableSoilwater.Table_Psd)
			end
			
			if option.psd.Table_Psd_θΨ_θ
				table.psd.PSD_θΨ_θ(IdSelect, N_iZ, hydroPsd, param, path.tableSoilwater.Table_Psd_θΨ_θ)
			end
		end # option.run.IntergranularMixingPsd

		if option.run.InfiltBest # <>=<>=<>=<>=<>
			table.infilt.HYDRO_INFILT(hydroInfilt, IdSelect, KunsatModel_Infilt, N_iZ, path.tableSoilwater.Table_HydroInfilt)

			table.infilt.INFILT(IdSelect, N_iZ, infiltOutput, path.tableSoilwater.Table_Infilt)
		end # option.run.InfiltBest

	# PRINT OUTPUT ======================================================================================
	# if option.other.Ploting && !option.run.Hypix
	# println("		=== START: PLOTTING  ===")
	
	# 	# if option.smap.Plot_Kunsat  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# 	# 	plotSmap.PLOT_KUNSAT(hydro, N_iZ, smap; N_Se= 1000)
	# 	# end

	# 	if option.run.HydroLabθΨ ≠ :No && option.hydro.Plot_θΨ # <>=<>=<>=<>=<>

	# 		if option.dataFrom.Smap
	# 			plotSmap.makie.HYDROPARAM(Ψ_θΨobs, Ψ_KΨobs, θ_θΨobs, N_θΨobs, N_iZ, N_KΨobs, K_KΨobs, IdSelect, hydro, Kₛ_Model, path; smap=smap)

	# 			# plotSmap.HYDROPARAM(Ψ_θΨobs, Ψ_KΨobs, θ_θΨobs, N_θΨobs, N_iZ, N_KΨobs, K_KΨobs, IdSelect, hydro, Kₛ_Model; N_Se=1000, smap=[])
	# 		else
	# 			plot.lab.HYDROPARAM(Ψ_θΨobs, Ψ_KΨobs, θ_θΨobs, N_θΨobs, N_iZ, N_KΨobs, K_KΨobs, IdSelect, hydro, Kₛ_Model, path.plotSoilwater.Plot_θΨK, path.option.ModelName)
	# 		end	
	# 	end # option.run.HydroLabθΨ
	# 	if option.run.IntergranularMixingPsd && option.psd.Plot_θr # <>=<>=<>=<>=<>
	# 		plot.psd.PLOT_θr(∑Psd, N_iZ, hydro, hydroPsd, path.plotSoilwater.Plot_Psd_θr, path.plotSoilwater.Plot_IMP_model)
	# 	end
	# 	if option.run.IntergranularMixingPsd && option.psd.Plot_IMP_Model # <>=<>=<>=<>=<>
	# 		plot.psd.PLOT_IMP_MODEL(IdSelect, Rpart, N_Psd, ∑Psd, Psd, N_iZ, hydro, paramPsd) 
	# 	end
	# 	if  option.run.IntergranularMixingPsd && option.psd.Plot_Psd_θΨ && !option.psd.HydroParam
	# 		println("			~ PSD WARNING Sorry cannot plot Plot_Psd_θΨ as option.psd.HydroParam==false ~")
	# 	end
	# 	if option.run.IntergranularMixingPsd && option.psd.Plot_Psd_θΨ && option.psd.HydroParam # <>=<>=<>=<>=<>
	# 		plot.psd.PLOT_PSD_θΨ(Ψ_θΨobs, Ψ_Rpart, θ_θΨobs, θ_Rpart, N_θΨobs, N_iZ, N_Psd, IdSelect, hydroPsd, hydro, path.plotSoilwater.Plot_Psd_θΨ)
	# 	end
	# 	if option.run.InfiltBest && option.infilt.Plot_∑Infilt  # <>=<>=<>=<>=<>
	# 		plot.infilt.PLOT_∑INFILT(IdSelect, N_Infilt, N_iZ, ∑Infilt_Obs, Tinfilt, ∑Infilt_3D, ∑Infilt_1D, infiltOutput, path.plotSoilwater.Plot_∑infilt_Opt )
	# 	end
	# 	# if option.run.InfiltBest && option.infilt.Plot_SeIni_Range # <>=<>=<>=<>=<>
	# 	# Removing GRUtils software to avoid conflict
	# 	# 	# plot.infilt.PLOT_∑INFILT_SEINI(hydroInfilt, IdSelect, infiltOutput, infiltParam, N_iZ)
	# 	# end

	# 	if option.run.InfiltBest && option.infilt.Plot_θΨ
	# 		if option.run.HydroLabθΨ ≠ :No
	# 			plot.infilt.PLOT_∑INFILT_θΨ(hydroInfilt, IdSelect, N_iZ, path.plotSoilwater.Plot_∑infilt_θΨ; hydro=hydro)
	# 		else
	# 			plot.infilt.PLOT_∑INFILT_θΨ(hydroInfilt, IdSelect, N_iZ, path.plotSoilwater.Plot_∑infilt_θΨ)
	# 		end # option.run.HydroLabθΨ
	# 	end # option.run.InfiltBest

	# println("=== END: PLOTTING  === \n")
	# end # if option.other.Ploting

	# Playing sounds...
		println("\007")

end  # function: START_TOOLBOX
# ..............................................................

println("\n\n ===== START SOIL WATER TOOLBOX =====")
	@time START_TOOLBOX()
println("==== END SOIL WATER TOOLBOX ====")