# =============================================================
#		MODULE: tool
# =============================================================
module tool
	# =============================================================
	#		module: normalize
	# =============================================================
	module norm
	
		export ∇NORM_2_PARAMETER

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : HYDRO_ADJEUSTMENTS
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∇NORM_2_PARAMETER(∇P, P_Min, P_Max)
				return P = ∇P * (P_Max - P_Min) + P_Min
			end  # function: HYDRO_ADJEUSTMENTS
		
	end  # module: normalize
	# ............................................................



	# =============================================================
	#		MODULE: readWrite
	# =============================================================
	module readWrite
		import DelimitedFiles
		export FIELDNAME_2_STRUCT_VECT, STRUCT_2_FIELDNAME, READ_HEADER, READ_ROW_SELECT, TOML_2_STRUCT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : READ_HEADER_FAST
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function READ_HEADER_FAST(Data, Header, Name)
			N_X, N_Y = size(Data) # Size of the array

			Header = reshape(Header, N_Y, 1)
			
		# Getting the column which matches the name of the header
			Name = replace(Name, " " => "") # Remove white spaces
			try
				iColumn = Int64(findfirst(isequal(Name), Header)[1])
				global Data_Output = Data[1:N_X,iColumn]
			catch
				println(Header)
				error("\n          SOILWATERTOOLBOX ERROR: cannot find  $Name   \n \n")
			end
		return Data_Output, N_X
		end # function READ_HEADER


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : READ_HEADER
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function READ_HEADER(Path, Name)
				# Read data
					Data =  DelimitedFiles.readdlm(Path, ',')
					N_X, N_Y = size(Data) # Size of the array
					
				# Reading header
					Header = fill("", N_Y)
					for i = 1:N_Y
						Header[i] = Data[1,i]
					end

				# Getting the column which matches the name of the header
					Name = replace(Name, " " => "") # Remove white spaces
	
					try
						global Data_Output = Data[2:N_X,findfirst(isequal(Name), Header)]
					catch
						println(Header)
						error("\n \n SOILWATERTOOLBOX ERROR: cannot find  $Name  in $Path \n \n")
					end

					N_X -= 1 # To take consideration of the header

			return Data_Output, N_X
			end # function READ_HEADER


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : READ_ROW_SELECT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function READ_ROW_SELECT(IdSelect::Vector{Int64}, Data, Header, Name::String, N_iZ::Int64; N_Point_Max=1000)

				Data_Output, N_X = READ_HEADER_FAST(Data, Header, Name)

			# ===========================================
			# Only keeping data which is selected
			# ===========================================
				Id_Data = Int64.(Data[1:end,1])

            N_Point = fill(0::Int64, N_iZ)
            iSelect = 1; iPoint = 1

				if isempty(Data_Output[1])
					Flag_String = false
					Data_Select = zeros(Union{Float64,Missing}, (N_iZ, N_Point_Max))

				elseif typeof(Data_Output[1]) == SubString{String}
					Data_Select = fill("A"::String, (N_iZ, N_Point_Max))
					Flag_String = true

				else
					Flag_String = false
					Data_Select = fill(0.0::Float64, (N_iZ, N_Point_Max))
				end

				# For all soils in the file
				for i = 1:N_X
					if Id_Data[i] == IdSelect[iSelect] # Only append Ids which correspond to the selected one
						if Flag_String
							Data_Output[i] = replace(Data_Output[i], " " => "") # Remove white spaces
						end
						if isempty(Data_Output[i])
							Data_Select[iSelect,iPoint] = missing
						else
							Data_Select[iSelect,iPoint] = Data_Output[i]
						end
						# append!(Data_Select, Data_Output[i])
						N_Point[iSelect] += 1
						iPoint += 1
					end

					# Since there are many Data_Output with the same Id only update IdSelect if we are changing soils and IdSelect[iSelect] == Id_Data[i]
					if i ≤ N_X -1
						if Id_Data[i+1] > Id_Data[i] && IdSelect[iSelect] == Id_Data[i] && iSelect ≤ N_iZ -1
							iSelect += 1
							iPoint = 1
						end # if:
					end # if: i ≤ N_X -1
				end # for: i = 1:N_X
	
		return Data_Select, N_Point
		end # function READ_ROW_SELECT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : FIELDNAME_2_STRUC
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function FIELDNAME_2_STRUCT_VECT(Structure, NameStruct)
				N_FieldName = length(fieldnames(Structure))

				FieldName_String = Array{Symbol}(undef, (N_FieldName))
				i = 1
				for FieldNames in fieldnames(Structure)
					FieldName_String[i] = FieldNames 
					i += 1
				end

				NameStruct.FieldName = FieldName_String

			return NameStruct
			end  # function: FIELDNAMES


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : STRUCT_2_FIELDNAMES
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function STRUCT_2_FIELDNAME(N_iZ, Structure)
				FieldName_Array = propertynames(Structure)

				N_FieldName = length(FieldName_Array)

				# Matrix
					Matrix = fill(0.0::Float64, (N_iZ, N_FieldName))

					i = 1
					for FieldName in FieldName_Array
						Struct_Array = getfield(Structure, FieldName)
		
						Matrix[1:N_iZ,i] .= Struct_Array
						i += 1
					end
				
				# HEADER
					FieldName_String = fill(""::String, N_FieldName)
					i=1
					for FieldNames in FieldName_Array
						FieldName_String[i] =  String(FieldNames)
						i += 1
					end
				return Matrix, FieldName_String
				end # function STRUCT_2_FIELDNAME

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TOML_2_STRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function TOML_2_STRUCT2(Structure, TomlParse)
      # LOOPING THROUGH THE DICT
      for (iKey, iValue₀) in TomlParse
      for iValue in (keys(iValue₀))
         if uppercase.(iKey) == (string(typeof(Structure)))
            setfield!(Structure, Symbol(iValue), TomlParse[iKey][iValue])
         end 
      end
   end
   return Structure
   end  # function: TOML_2_STRUCT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TOML_2_STRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TOML_2_STRUCT(Structure, TomlParse; MyType_LowerCase=true, MyType=:MyType)
			if MyType_LowerCase == false
				MyType = string(MyType)
			else
				MyType = lowercase.(string(Structure))
			end

			Output = NamedTuple{Tuple(Symbol.(keys(TomlParse[MyType])))}(values(TomlParse[MyType]))
		return Structure(Output...)
		end # function TOML_2_STRUC


	end  # module readWrite ************************
	# ............................................................
end  # module tool
# ............................................................