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
		import CSV, Tables
		export FIELDNAME_2_STRUCT_VECT, STRUCT_2_FIELDNAME, READ_HEADER, READ_ROW_SELECT, READ_θΨK_2D, READ_STRUCT_SIMPLE

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : READ_HEADER_FAST
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function READ_HEADER_FAST(Data, Header, Name)
				N_X, N_Y = size(Data) # Size of the array
				
				# Getting the column which matches the name of the header
				Name = replace(Name, " " => "") # Remove white spaces
				
				# Data_Output = []
				# try
					Header = reshape(Header, N_Y, 1)
					iColumn = Int64(findfirst(isequal(Name), reshape(Header, N_Y, 1))[1])

					Data_Output =  Data[1:N_X,iColumn]

					# convert(Vector{any}, Data_Output)

				# catch
				# 	println(Header)
				# 	error("\n          SOILWATERTOOLBOX ERROR: cannot find  $Name   \n \n")
				# end
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
			function READ_ROW_SELECT(IdSelect::Vector{Int64}, Data, Header, Name::String, NiZ::Int64; N_Point_Max=1000)

				Data_Output, N_X = READ_HEADER_FAST(Data, Header, Name)

				# if NiZ ≠ N_X
				# 	error("READ_ROW_SELECT  NiZ ≠ N_X")
				# end

			# ===========================================
			# Only keeping data which is selected
			# ===========================================
				Id_Data = Int64.(Data[1:end,1])

				# Determening if there is existence of an empty cell
				Flag_IsEmpty = false
				for i = 1:N_X
					if isempty(Data_Output[i])
						Flag_IsEmpty = true
						break
					end
				end

				if Flag_IsEmpty
					Flag_String = false
					Data_Select = zeros(Union{Float64,Missing}, (NiZ, N_Point_Max))

				elseif typeof(Data_Output[1]) == SubString{String}
					Data_Select = fill(Data_Output[1], (NiZ, N_Point_Max))
					Flag_String = true

				else
					Flag_String = false
					Data_Select = fill(0.0::Float64, (NiZ, N_Point_Max))
				end

				# For all soils in the file
				iSelect = 1; iPoint = 1
            N_Point = fill(0::Int64, NiZ)
				
				for i = 1:N_X
					if Id_Data[i] > IdSelect[iSelect] && iSelect ≠ NiZ
						error("READ_ROW_SELECT problem with no matching id:  i= $i IdSelect[iSelect] = $( IdSelect[iSelect])< Id_Data[i] = $(Id_Data[i])")
					end

					if Id_Data[i] == IdSelect[iSelect] # Only append Ids which correspond to the selected one
						if Data_Output[i] == SubString{String}
							Data_Output[i] = replace(Data_Output[i], " " => "") # Remove white spaces
						end
						if isempty(Data_Output[i])
							Data_Select[iSelect,iPoint] = missing
						else
							Data_Select[iSelect,iPoint] = Data_Output[i]
						end
						N_Point[iSelect] += 1
		
						if (i ≤ N_X - 1) && (iSelect ≤ NiZ -1) && (Id_Data[i+1] > Id_Data[i]) 
							iSelect += 1
							iPoint = 1
						else
							iPoint += 1
						end # if:
					end # if: i ≤ N_X -1
				end

				if iSelect ≠ NiZ
					error("READ_ROW_SELECT error: iSelect=$iSelect ≠ NiZ=$NiZ")
				end

					# Since there are many Data_Output with the same Id only update IdSelect if we are changing soils and IdSelect[iSelect] == Id_Data[i]
					# if i ≤ N_X - 1
					# 	if (Id_Data[i+1] > Id_Data[i]) && (IdSelect[iSelect] == Id_Data[i]) && (iSelect ≤ NiZ -1)
					# 		iSelect += 1
					# 		iPoint = 1
					# 	end # if:
					# end # if: i ≤ N_X -1
				# end # for: i = 1:N_X
	
		return Data_Select, N_Point
		end # function READ_ROW_SELECT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : READ_θΨK_2D
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function READ_θΨK_2D(Data, Header, IdSelect, NiZ)
				Ψheader = fill(0.0, length(Header)-1)
				Xdata = zeros(Union{Float64,Missing}, (NiZ, length(Header)-1))

				iΨ=1
				for iHeader_Ψ in Header
					if iHeader_Ψ  ≠ "\ufeffId" && iHeader_Ψ  ≠ "Id"
						Xdata_Column, ~ = readWrite.READ_ROW_SELECT(IdSelect, Data, Header, string(iHeader_Ψ), NiZ; N_Point_Max=1)

						# Reading Xdata for every column
							Xdata[1:NiZ, iΨ] = Xdata_Column[1:NiZ]

						# Converting the value of header into number by removing all characters Ψ (e.g. H_100_mm) -> 100
							iFindall = findall("_", iHeader_Ψ)
							iHeader_Ψ =  iHeader_Ψ[iFindall[1][1]+1:iFindall[2][1]-1]
							iHeader_Float =  parse(Float64, iHeader_Ψ)
							Ψheader[iΨ] = iHeader_Float
							iΨ += 1
					end # if iHeader_Ψ  ≠ "\ufeffId" && iHeader_Ψ  ≠ "Id"
				end # for iHeader_Ψ in Header

				Ψdata_NoMissing = fill(0.0::Float64, (NiZ, length(Header)-1))
				Xdata_NoMissing = fill(0.0::Float64, (NiZ, length(Header)-1))
				NΨobs = fill(0::Int64, NiZ)

				for iZ=1:NiZ
					Header_Ψ₂ = deepcopy(Ψheader)
					
					# Finding Missing data
						Ψheader_Missing = findall(x->x===missing, Xdata[iZ, 1:length(Header)-1])

					# Deleating the missing data in a row
						Xdata_Row = deleteat!(Xdata[iZ, 1:length(Header)-1], Ψheader_Missing)

						deleteat!(Header_Ψ₂, Ψheader_Missing)

					# Finding the number of data points for every iZ which are not missing
						NΨobs[iZ] = length(Header_Ψ₂)

					# Writing non missing data
						Xdata_NoMissing[iZ, 1:NΨobs[iZ]] = Xdata_Row[1:NΨobs[iZ]]

						Ψdata_NoMissing[iZ,1:NΨobs[iZ]] = Header_Ψ₂[1:NΨobs[iZ]]

				end # for iZ=1:NiZ	
			return NΨobs, Xdata_NoMissing, Ψdata_NoMissing
			end  # function: Read_2	


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
			function STRUCT_2_FIELDNAME(N::Int64, Structure)
				FieldName_Array = propertynames(Structure)

				N_FieldName = length(FieldName_Array)

				# Matrix
					Matrix = fill(0.0, (N, N_FieldName))

					i = 1
					for FieldName in FieldName_Array
						Struct_Array = getfield(Structure, FieldName)
						if isa(Struct_Array, Array)
							Matrix[1:N,i] = Float64.(Struct_Array[1:N])
						else
							Matrix[1,i] = Float64.(Struct_Array)
						end
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
		#		FUNCTION : READ_STRUCT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function READ_STRUCT_SIMPLE(structure, Path::String)
				# println("    ~  $(Path) ~")

				# Read data
					Data = CSV.File(Path, header=true)
				
				# Select data of interest
					NiZ = size(Data)[1] # Initial
			
				# Reading the Model data
				for iFieldname in propertynames(structure)			
					
					Output_Vector = convert(Vector{Float64}, Tables.getcolumn(Data, iFieldname))

					if typeof(getfield(structure, iFieldname)) == Vector{Float64}
						setfield!(structure, Symbol(iFieldname), Float64.(Output_Vector))
					else
						setfield!(structure, Symbol(iFieldname), Float64(Output_Vector[1]))
					end
				end

			return structure, NiZ
			end  # function: READ_STRUCT
		#----------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : READ_STRUCT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function READ_STRUCT2(structure, Path::String; iStart=1, iEnd=2^63 - 1)
				println("    ~  $(Path) ~")

				# Read data
					Data, Header = DelimitedFiles.readdlm(Path, ',', header=true, use_mmap=true)

				# Select data of interest
					NiZ = size(Data)[1] # Initial
					iEnd= min(NiZ, iEnd)
					Data = Data[iStart:iEnd,1:end]
					NiZ = iEnd - iStart + 1 # Final

				# Reading the Model data
				for iFieldname in propertynames(structure)

					# Putting the values of Output_Vector into structure
					Output_Vector = fill(0.0::Float64, NiZ)					
					try
						Output_Vector, Ndata = readWrite.READ_HEADER_FAST(Data, Header, string(iFieldname))
					catch
						# @warn "SoilWater-ToolBox: cannong find $iFieldname"
						Output_Vector = fill(0.0::Float64, NiZ)
					end

					try
						setfield!(structure, Symbol(iFieldname), Float64.(Output_Vector))
					catch
						setfield!(structure, Symbol(iFieldname), Float64(Output_Vector[1]))
					end
				end

			return structure, NiZ
			end  # function: READ_STRUCT
	#----------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TOML_2_STRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # function TOML_2_STRUCT2(Structure, TomlParse)
   #    # LOOPING THROUGH THE DICT
   #    for (iKey, iValue₀) in TomlParse
   #    for iValue in (keys(iValue₀))
   #       if uppercase.(iKey) == (string(typeof(Structure)))
   #          setfield!(Structure, Symbol(iValue), TomlParse[iKey][iValue])
   #       end 
   #    end
   # end
   # return Structure
   # end  # function: TOML_2_STRUCT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TOML_2_STRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function TOML_2_STRUCT(Structure, TomlParse; MyType_LowerCase=true, MyType=:MyType)
		# 	if MyType_LowerCase == false
		# 		MyType = string(MyType)
		# 	else
		# 		MyType = lowercase.(string(Structure))
		# 	end

		# 	Output = NamedTuple{Tuple(Symbol.(keys(TomlParse[MyType])))}(values(TomlParse[MyType]))
		# return Structure(Output...)
		# end # function TOML_2_STRUC


	end  # module readWrite ************************
	# ............................................................
end  # module tool
# ............................................................