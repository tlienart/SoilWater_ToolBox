	# NOT READY YET
	
	
	# _______________________ START: ChangeHydroModel _______________________ 
	if option.run.ChangeHydroModel
		# STRUCTURES
			NiZ=1000
			hydroTranslate = hydroStruct.HYDROSTRUCT(option.hydro, NiZ)
			hydroOtherTranslate = hydroStruct.HYDRO_OTHERS(NiZ)
			hydroTranslate, optimTranslate = reading.HYDRO_PARAM(option.hydro, hydroTranslate, NiZ, path.inputGuiSoilwater.GUI_HydroParam)
		
		# Temporary Id
			IdSelect = collect(1:1:NiZ)
	
		# Deriving a table of θ(Ψ)
			table.hydroLab.TABLE_EXTRAPOINTS_θΨ(hydroTranslate, IdSelect, NiZ, path.inputSoilwater.Ψθ, param.hydro.TableComplete_θΨ; Orientation="Vertical")
		
		# Deriving a table of K(θ)
			table.hydroLab.TABLE_EXTRAPOINTS_Kθ(hydroTranslate, IdSelect, param.hydro.TableComplete_θΨ, hydroTranslate.Ks[1:NiZ], NiZ::Int64, path.inputSoilwater.Kunsat)

		# Creating an Id output required by the program
			table.TABLE_ID(NiZ::Int64, path, path.inputSoilwater.IdSelect)
	end # Option
	# ------------------------END: ChangeHydroModel---------------------------  	