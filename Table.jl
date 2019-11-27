# =============================================================
#		MODULE: table
# =============================================================
module table

	# =============================================================
	#		MODULE: name
	# =============================================================
	module hydroParam
		import ...path, ...tool
		import DelimitedFiles
		export θΨK

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θΨK(Id_Select, N_SoilSelect, hydro)

			Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect, hydro)
			
			pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
			Matrix =  round.(Matrix, digits=3)
			open(path.Table_θΨK, "w") do io
				DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
				DelimitedFiles.writedlm(io, [Id_Select Matrix], ",")
			end
	 
			return
		end  # function:  θΨK
		
	end  # module hydro
	# ............................................................


	# =============================================================
	#		MODULE: psd
	# =============================================================
	module psd
		import ...path, ...tool
		import DelimitedFiles
		export PSD, θΨK_PSD

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PSD
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PSD(Id_Select, N_SoilSelect, paramPsd)
				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect,  paramPsd)
				
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

				open(path.Table_Psd, "w") do io
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [Int64.(Id_Select) round.(Matrix,digits=3)], ",")
				end
				return
			end


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK_PSD
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θΨK_PSD(Id_Select, N_SoilSelect, hydroPsd)
				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect, hydroPsd)
				
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
				Matrix =  round.(Matrix, digits=3)
				open(path.Table_θΨK_Psd, "w") do io
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [Id_Select Matrix], ",")
				end
				return
			end  # function:  θΨK
		
	end  # module psd
	# ............................................................
	
end  # module table
# ............................................................