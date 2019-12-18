# =============================================================
#		MODULE: read
# =============================================================
module read
	import ..option, ..path, ..cst, ..tool
	export ID, θΨ, KUNSATΨ, INFILTRATION, PSD, READ_ROW_SELECT, ρ_Ψθ, ρ_INFILTRATION, ρ_Psd

	mutable struct INFILT
		RingRadius
		θ_Ini
		γ
		β
	end


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ID
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ID()
			Id_True, N_iSoil_All = tool.readWrite.READ_HEADER(path.Id_Select, "TRUE")
			N_SoilSelect = sum(Id_True)

			Id_Select = Array{Int64}(undef, N_SoilSelect)
			iTrue = 1
			for iSoil in 1:N_iSoil_All
				if Id_True[iSoil] == 1
					Id_Select[iTrue] = iSoil
					iTrue += 1
				end	# Id_Soil == 1	
			end  # for: Id_Soil in Id_True
			
			return Id_Select, N_SoilSelect
		end  # function: ID


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θΨ(Id_Select, N_SoilSelect)
			Ψ_θΨ, N_θΨ = tool.readWrite.READ_ROW_SELECT(path.Ψθ, "H[mm]", Id_Select, N_SoilSelect)
			
         θ_θΨ, ~    = tool.readWrite.READ_ROW_SELECT(path.Ψθ, "Theta", Id_Select, N_SoilSelect)
			return θ_θΨ, Ψ_θΨ, N_θΨ
		end  # function: θΨ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ρ_Ψθ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρ_Ψθ(Id_Select, N_SoilSelect)
			
			ρbSoil_θΨ, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "ρbSoil[g_cm-3]", Id_Select, N_SoilSelect)
			
			ρp_Fine_θΨ, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "ρp_Fine[g_cm-3]", Id_Select, N_SoilSelect)

			ρ_Rock_θΨ, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "ρ_Rock[g_cm-3]", Id_Select, N_SoilSelect)
			
         RockW_θΨ, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "Rock%[g_g-3]", Id_Select, N_SoilSelect)

			return RockW_θΨ, ρ_Rock_θΨ, ρbSoil_θΨ, ρp_Fine_θΨ
		end # function: ρ_Ψθ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KUNSATΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KUNSATΨ(Id_Select, N_SoilSelect)
			Ψ_KΨ, N_KΨ = tool.readWrite.READ_ROW_SELECT(path.Kunsat, "H[mm]", Id_Select, N_SoilSelect)
			
         K_KΨ, ~    = tool.readWrite.READ_ROW_SELECT(path.Kunsat, "Kunsat[mm_s]", Id_Select, N_SoilSelect)
			return K_KΨ, Ψ_KΨ, N_KΨ 
		end  # function: θΨ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION(Id_Select, N_SoilSelect)
			Tinfilt, N_Infilt = tool.readWrite.READ_ROW_SELECT(path.Infiltration, "Tinfilt[s]", Id_Select, N_SoilSelect)
			
         ∑Infilt_Obs , ~   = tool.readWrite.READ_ROW_SELECT(path.Infiltration, "Cumul_Infiltration[mm]", Id_Select, N_SoilSelect)

			RingRadius , ~    = tool.readWrite.READ_ROW_SELECT(path.Infiltration_Param, "RingRadius[mm]", Id_Select, N_SoilSelect, N_Point_Max=1)
			
			θ_Ini , ~         = tool.readWrite.READ_ROW_SELECT(path.Infiltration_Param, "Theta_Ini[-]", Id_Select, N_SoilSelect, N_Point_Max=1)
			
			γ , ~             = tool.readWrite.READ_ROW_SELECT(path.Infiltration_Param, "Lambda[-]", Id_Select, N_SoilSelect, N_Point_Max=1)
			
         β , ~             = tool.readWrite.READ_ROW_SELECT(path.Infiltration_Param, "Beta[-]", Id_Select, N_SoilSelect, N_Point_Max=1)

			infiltParam = INFILT(RingRadius, θ_Ini, γ, β)

			return Tinfilt, ∑Infilt_Obs, N_Infilt, infiltParam
		end  # function: INFILTRATION


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ρ_PSD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρ_INFILTRATION(Id_Select, N_SoilSelect)

			ρbSoil_Infilt, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "ρbSoil[g_cm-3]", Id_Select, N_SoilSelect)
			
			ρp_Fine_Infilt, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "ρp_Fine[g_cm-3]", Id_Select, N_SoilSelect)

			ρ_Rock_Infilt, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "ρ_Rock[g_cm-3]", Id_Select, N_SoilSelect)
			
         RockW_Infilt, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "Rock%[g_g-3]", Id_Select, N_SoilSelect)

			return RockW_Infilt, ρ_Rock_Infilt, ρbSoil_Infilt, ρp_Fine_Infilt 
		end # function: ρ_PSD


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PSD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PSD(Id_Select, N_SoilSelect) # TODO make sure that the particles are ordered from smalest to largest
			Diameter_Psd, N_Psd 	= tool.readWrite.READ_ROW_SELECT(path.Psd, "Diameter[mm]", Id_Select, N_SoilSelect)
			∑Psd , ~ 				= tool.readWrite.READ_ROW_SELECT(path.Psd, "Cumul_Psd", Id_Select, N_SoilSelect)

			Rpart = @. Diameter_Psd / 2.0

			return Rpart, ∑Psd, N_Psd
		end  # function: PSD


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ρ_PSD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρ_PSD(Id_Select, N_SoilSelect)

			ρbSoil_Psd, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "ρbSoil[g_cm-3]", Id_Select, N_SoilSelect)
			
			ρp_Fine_Psd, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "ρp_Fine[g_cm-3]", Id_Select, N_SoilSelect)

			ρ_Rock_Psd, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "ρ_Rock[g_cm-3]", Id_Select, N_SoilSelect)
			
         RockW_Psd, ~ = tool.readWrite.READ_ROW_SELECT(path.ρ_Ψθ, "Rock%[g_g-3]", Id_Select, N_SoilSelect)

			return RockW_Psd, ρ_Rock_Psd, ρbSoil_Psd, ρp_Fine_Psd
		end # function: ρ_PSD
	
end  # module: read
# ............................................................		