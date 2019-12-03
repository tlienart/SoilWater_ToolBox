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
	export FIELDNAME_2_STRUCT, STRUCT_2_FIELDNAME
		
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