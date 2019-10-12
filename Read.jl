# =============================================================
#		MODULE: read
# =============================================================
module read
	using ..option, ..path, ..cst
	import DelimitedFiles

	export ID, KUNSATΨ, INFILTRATION, PSD, READ_ROW_SELECT

	mutable struct INFILT
		RingRadius
		Theta_Ini
		θs
		γ
		β
	end


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ID
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ID()
			Id_True, N_iSoil_All = READ_HEADER(path.Id_Select, "TRUE")
			N_SoilSelect = sum(Id_True)

			Id_Select = Array{Int64}(undef, N_SoilSelect)
			iTrue = 1
			for iSoil in 1:N_iSoil_All
				if Id_True[iSoil] == 1
					Id_Select[iTrue] = iSoil
					iTrue += 1
				end	# Id_Soil == 1	
			end  # for: Id_Soil in Id_True
			
			return Id_Select, Id_True, N_SoilSelect
		end  # function: ID


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θΨ(Id_True, N_SoilSelect)
			Ψ_θΨ, N_θΨ 	= READ_ROW_SELECT(path.Ψθ, "H[kPa]", Id_True, N_SoilSelect)
			Ψ_θΨ = Ψ_θΨ * cst.kPa_2_mm
			θ_θΨ, ~ 	= READ_ROW_SELECT(path.Ψθ, "Theta", Id_True, N_SoilSelect)
			return θ_θΨ, Ψ_θΨ, N_θΨ
		end  # function: θΨ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KUNSATΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KUNSATΨ(Id_True, N_SoilSelect)
			Ψ_KΨ, N_KΨ 	= READ_ROW_SELECT(path.Kunsat, "H[kPa]", Id_True, N_SoilSelect)
			K_KΨ, ~ 	= READ_ROW_SELECT(path.Kunsat, "Kunsat[mm_s]", Id_True, N_SoilSelect)
			return K_KΨ, Ψ_KΨ, N_KΨ 
		end  # function: θΨ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION(Id_True, N_SoilSelect)
			T, N_Infilt 	= READ_ROW_SELECT(path.Infiltration, "T[s]", Id_True, N_SoilSelect)
			∑Infilt , ~ 	= READ_ROW_SELECT(path.Infiltration, "Cumul_Infiltration[mm]", Id_True, N_SoilSelect)

			RingRadius , ~ 	=  READ_ROW_SELECT(path.Infiltration_Param, "RingRadius[mm]", Id_True, N_SoilSelect, N_Point_Max=1)
			Theta_Ini , ~ 	=  READ_ROW_SELECT(path.Infiltration_Param, "Theta_Ini[-]", Id_True, N_SoilSelect, N_Point_Max=1)
			θs, ~ 			=  READ_ROW_SELECT(path.Infiltration_Param, "Thetas_Bulk[-]", Id_True, N_SoilSelect, N_Point_Max=1)
			γ , ~ 			=  READ_ROW_SELECT(path.Infiltration_Param, "Lambda[-]", Id_True, N_SoilSelect, N_Point_Max=1)
			β , ~ 			=  READ_ROW_SELECT(path.Infiltration_Param, "Beta[-]", Id_True, N_SoilSelect, N_Point_Max=1)

			infilt = INFILT(RingRadius, Theta_Ini, θs, γ, β)

			return T, ∑Infilt, N_Infilt, infilt
		end  # function: INFILTRATION


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PSD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PSD(Id_True, N_SoilSelect)
			Diameter, N_Psd 	= READ_ROW_SELECT(path.Psd, "Diameter[mm]", Id_True, N_SoilSelect)
			∑Psd , ~ 	= READ_ROW_SELECT(path.Psd, "Cumul_Psd", Id_True, N_SoilSelect)
			return Diameter, ∑Psd, N_Psd
		end  # function: PSD


	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : READ_HEADER
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function READ_HEADER(Path, Name)
			# Read data
				Data =  DelimitedFiles.readdlm(Path, ',')
				N_X, N_Y = size(Data) # Size of the array

				# Reading header
				Header = fill("", N_Y)
				for i in 1:N_Y
					Header[i] = Data[1,i]
				end

				# Getting the column which matches the name of the header
				Name = replace(Name, " " => "") # Remove white spaces
				try
					global Data_Output = Data[2:N_X,findfirst(isequal(Name), Header)]
				catch
					error("\n \n SOILWATERTOOLBOX ERROR: cannot find   $Name  in $Path \n \n")
				end

				N_X -= 1 # To take consideration of the header

				return Data_Output, N_X
			end # function READ_HEADER


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : READ_HEADER_VERTICAL_FLAG
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function READ_ROW_SELECT(Path::String, Name::String, Id_Select::Vector{Int64}, N_SoilSelect::Int64; N_Point_Max=500)
			# READ DATA
				Data =  DelimitedFiles.readdlm(Path, ',')
				N_X, N_Y = size(Data) # Size of the array

				# Get the ID of the data
				Id_Data = Int64.(Data[2:N_X,1])

			# READ HEADER
				Header = fill("", N_Y)
				for i in 1:N_Y
					Header[i] = Data[1,i]
				end

				# Getting the column which matches the name of the header
				Name = replace(Name, " " => "") # Remove white spaces
		
				Data_Output = Float64.(Data[2:N_X,findfirst(isequal(Name), Header)])

				N_X -= 1 # Take consideration of the header

			# ===========================================
			# Only keeping data which is selected
			# ===========================================
				# Data_Select = zeros(Float64, (N_SoilSelect, N_Point_Max) )
				Data_Select = Array{Float64}(undef, (N_SoilSelect, N_Point_Max))
				N_Point = zeros(Int64, N_SoilSelect) 
				
				iSelect = 1; iPoint = 1
				
				# For all soils in the file
				for i in 1:N_X
					if Id_Data[i] == Id_Select[iSelect] # Only append Ids which correspond to the selected one
						Data_Select[iSelect,iPoint] =  Data_Output[i]
						# append!(Data_Select, Data_Output[i])
						N_Point[iSelect] += 1
						iPoint += 1
					end

					# Since there are many Data_Output with the same Id only update Id_Select if we are changing soils and Id_Select[iSelect] == Id_Data[i]
					if i ≤ N_X -1
						if Id_Data[i+1] > Id_Data[i] && Id_Select[iSelect] == Id_Data[i]
							iSelect += 1
							iPoint = 1
						end # if:
					end # if: i ≤ N_X
				end # for: i in 1:N_X

				return Data_Select, N_Point
		end # function READ_ROW_SELECT

end  # module: read
# ............................................................		