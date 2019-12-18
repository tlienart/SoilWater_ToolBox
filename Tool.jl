# =============================================================
#		MODULE: tool
# =============================================================
module tool

	# =============================================================
	#		MODULE: array
	# =============================================================
	module array
		function SEARCH_INDEX(Array, SearchValue)
			N = length(Array)
			iSearchValue = 1
			Value_SearchValue=1.
			Err_2 = 100000000000000.
			
			for i in 1:N
				Err_1 = abs(Array[i] - SearchValue)
			
				if Err_1 < Err_2
					iSearchValue = i
					Value_SearchValue = Array[i]
					Err_2 = Err_1
				end
			end
			return iSearchValue
		end # function SEARCH_INDEX
	end  # module: array
	# ............................................................


	# =============================================================
	#		MODULE: readWrite
	# =============================================================
	module readWrite
		import DelimitedFiles
		export FIELDNAME_2_STRUCT, STRUCT_2_FIELDNAME

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
			function READ_ROW_SELECT(Path::String, Name::String, Id_Select::Vector{Int64}, N_SoilSelect::Int64; N_Point_Max=50)
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
							if Id_Data[i+1] > Id_Data[i] && Id_Select[iSelect] == Id_Data[i] && iSelect ≤ N_SoilSelect -1
								iSelect += 1
								iPoint = 1
							end # if:
						end # if: i ≤ N_X
					end # for: i in 1:N_X

					return Data_Select, N_Point
			end # function READ_ROW_SELECT

		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : FIELDNAME_2_STRUC
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function FIELDNAME_2_STRUCT(Structure, NameStruct)
				N_FieldName = length(fieldnames(Structure))

				FieldName_String = Array{Symbol}(undef, (N_FieldName))
				i = 1
				for FieldNames in fieldnames(Structure)
					FieldName_String[i] = FieldNames 
					i += 1
				end

				NameStruct.FieldName = FieldName_String

				# for i=1:N_FieldName
				# 	Struct_Array = getfield(NameStruct, FieldName_String[i])
				# end

				return NameStruct
			end  # function: FIELDNAMES


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : STRUCT_2_FIELDNAMES
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function STRUCT_2_FIELDNAME(N_SoilSelect, Structure)
				N_FieldName = length(Structure.FieldName) - 1

				Matrix = Array{Float64}(undef, (N_SoilSelect, N_FieldName))
				
				for i=1:N_FieldName
					Struct_Array = getfield(Structure, Structure.FieldName[i])
					Matrix[1:N_SoilSelect,i] = Struct_Array[1:N_SoilSelect]
				end

				FieldName_String = Array{String}(undef, N_FieldName)
				i=1
				for FieldNames in Structure.FieldName
					FieldName_String[i] =  String(FieldNames)
					if i == N_FieldName
						break
					end
					i += 1
				end

				return Matrix, FieldName_String
			end # function STRUCT_2_FIELDNAME
		
	end  # module readWrite
	# ............................................................
	
end  # module tool
# ............................................................