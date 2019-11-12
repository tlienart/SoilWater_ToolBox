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
		function θΨK(Id_Select, Of, Of_θΨ, Of_Kunsat, N_SoilSelect, hydro)

			Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect, hydro)
			
			pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

			Nse = 1 .- Of
			Nse_θΨ = 1 .- Of_θΨ
			Nse_Kunsat = 1 .- Of_Kunsat

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
		export PSD, PSD_θr

		function PSD(Id_Select, N_SoilSelect, psdparam)
			Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect,  psdparam)
			
			pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

			open(path.Table_Psd, "w") do io
				DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
				DelimitedFiles.writedlm(io, [Int64.(Id_Select) round.(Matrix,digits=3)], ",")
			end
		end


		function PSD_θr(Err_θr_Psd, Id_Select, N_SoilSelect, θr, θr_Psd)
			open(path.Table_θr, "w") do io
				Header = ["Id", "θr_Obs", "θr_Psd", "Err_θr_Psd"]
				DelimitedFiles.writedlm(io,[Header] , ",",) # Header

				DelimitedFiles.writedlm(io, [Int64.(Id_Select) round.(θr,digits=4) round.(θr_Psd,digits=4)  round.(Err_θr_Psd,digits=3)], ",")
			end
		end

		
	end  # module psd
	# ............................................................






	

		




		





	
end  # module table
# ............................................................