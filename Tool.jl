# =============================================================
#		MODULE: tool
# =============================================================
module tool

	export STRUCT_2_FIELDNAME
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : STRUCT_2_FIELDNAMES
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function STRUCT_2_FIELDNAME(Structure, NameStruct)
			N_FieldName = length(fieldnames(Structure))

			FieldName_String = Array{Symbol}(undef, (N_FieldName))
			i = 1
			for FieldNames in fieldnames(Structure)
				FieldName_String[i] = FieldNames 
				i += 1
			end

			NameStruct.FieldName = FieldName_String

			for i=1:N_FieldName
				Struct_Array = getfield(NameStruct, FieldName_String[i])
			end

			return NameStruct
		end  # function: FIELDNAMES
end  # module tool
# ............................................................